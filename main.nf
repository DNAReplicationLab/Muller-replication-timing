#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * Muller-replication-timing
 *
 * DNA copy-number measurement of genome replication dynamics by high-throughput
 * sequencing: the sort-seq, sync-seq and MFA-seq family.
 *
 * Reimplementation of the pipeline described in
 *   Batrakou, Muller, Wilson & Nieduszynski (2020) Nature Protocols 15:1255-1284
 *   doi:10.1038/s41596-019-0287-7
 * which combined DNAReplicationLab/localMapper (mapping and binning) with
 * DNAReplicationLab/Repliscope (ratios, normalisation, Trep).
 */

include { SRA_FETCH          } from './modules/sra_fetch.nf'
include { BOWTIE2_INDEX      } from './modules/bowtie2_index.nf'
include { BOWTIE2_ALIGN      } from './modules/bowtie2_align.nf'
include { COUNT_BINS         } from './modules/count_bins.nf'
include { REPLISCOPE_QC      } from './modules/repliscope.nf'
include { REPLISCOPE_PROFILE } from './modules/repliscope.nf'
include { REPLISCOPE_TREP    } from './modules/repliscope.nf'
include { BED_TO_BIGWIG as BIGWIG_SMOOTH } from './modules/bigwig.nf'
include { BED_TO_BIGWIG as BIGWIG_TREP   } from './modules/bigwig.nf'

def helpMessage() {
    log.info """
    Muller-replication-timing  --  sort-seq / sync-seq / MFA-seq replication profiles

    Usage:
      nextflow run main.nf --samplesheet SHEET.csv --comparisons COMPS.csv \\
                           --fasta genome.fa --outdir results

    Required:
      --samplesheet   CSV: sample,run,platform,layout,variant,role,strain,timepoint,series
                      'run' is an SRA run accession, or leave blank and give --fastq_dir
      --comparisons   CSV: comparison,rep,nonRep,variant,norm_mode,norm_value,description
      --fasta         Reference genome FASTA (default: UCSC sacCer3)

    Key options:
      --bin_size            Bin width in bp (default ${params.bin_size})
      --bin_lower_limit     Drop bins below this fraction of the genome median (default ${params.bin_lower_limit})
      --fastq_dir           Use local FASTQ files instead of fetching from SRA
      --outdir              Results directory (default ${params.outdir})
    """.stripIndent()
}

workflow {
    if (params.help) {
        helpMessage()
        exit(0)
    }

    // ---- inputs ----------------------------------------------------------
    ch_samples = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            assert row.platform in ['ILLUMINA', 'ABI_SOLID'] : "unknown platform: ${row.platform}"
            assert row.layout   in ['single', 'paired']      : "unknown layout: ${row.layout}"
            assert row.role     in ['rep', 'nonRep']         : "unknown role: ${row.role}"
            row
        }

    if (params.bowtie2_index) {
        // An empty glob stages nothing and Bowtie 2 then fails deep inside a pipe with
        // an unhelpful samtools error, so check the index is really there up front.
        def idx_files = files("${params.bowtie2_index}*.bt2")
        if (idx_files.size() < 6) {
            error "--bowtie2_index ${params.bowtie2_index} matched ${idx_files.size()} .bt2 " +
                  "files; a Bowtie 2 index needs 6. Check the path and basename, or omit " +
                  "--bowtie2_index to build one from --fasta."
        }
        ch_index = Channel.value(tuple(idx_files, file(params.bowtie2_index).name))
    }
    else {
        ch_index = BOWTIE2_INDEX(Channel.fromPath(params.fasta, checkIfExists: true)).index
    }

    // ---- reads -----------------------------------------------------------
    if (params.fastq_dir) {
        ch_reads = ch_samples.map { meta ->
            def files = meta.layout == 'paired'
                ? [file("${params.fastq_dir}/${meta.sample}_1.fastq.gz", checkIfExists: true),
                   file("${params.fastq_dir}/${meta.sample}_2.fastq.gz", checkIfExists: true)]
                : [file("${params.fastq_dir}/${meta.sample}.fastq.gz", checkIfExists: true)]
            tuple(meta, files)
        }
    } else {
        ch_reads = SRA_FETCH(ch_samples).reads
    }

    // ---- map and count ---------------------------------------------------
    BOWTIE2_ALIGN(ch_reads, ch_index)
    COUNT_BINS(BOWTIE2_ALIGN.out.bam)
    REPLISCOPE_QC(COUNT_BINS.out.counts)

    // ---- pair up each comparison with its two count files ----------------
    ch_counts = COUNT_BINS.out.counts.map { meta, bed -> tuple(meta.sample, bed) }

    ch_comparisons = Channel
        .fromPath(params.comparisons, checkIfExists: true)
        .splitCsv(header: true)

    ch_profile_in = ch_comparisons
        .map { comp -> tuple(comp.rep, comp) }
        .combine(ch_counts, by: 0)
        .map { repName, comp, repBed -> tuple(comp.nonRep, comp, repBed) }
        .combine(ch_counts, by: 0)
        .map { nonRepName, comp, repBed, nonRepBed -> tuple(comp, repBed, nonRepBed) }

    REPLISCOPE_PROFILE(ch_profile_in)

    // ---- Trep from the sync-seq time course ------------------------------
    ch_sync = REPLISCOPE_PROFILE.out.ratio
        .filter { comp, tsv -> comp.variant == 'syncseq' }
        .map { comp, tsv -> tuple(comp.comparison, tsv) }
        .toSortedList { a, b ->
            // order by the timepoint embedded in the comparison name
            (a[0] =~ /(\d+)min/)[0][1].toInteger() <=> (b[0] =~ /(\d+)min/)[0][1].toInteger()
        }

    // Trep needs at least two timepoints to fit; skip the step entirely otherwise
    // (e.g. a run covering only sort-seq or MFA-seq samples).
    ch_sync_ok    = ch_sync.filter { list -> list.size() >= 2 }
    ch_trep_files = ch_sync_ok.map { list -> list.collect { it[1] } }
    ch_trep_times = ch_sync_ok.map { list -> list.collect { (it[0] =~ /(\d+)min/)[0][1].toInteger() } }

    REPLISCOPE_TREP(ch_trep_files, ch_trep_times)

    // ---- browser tracks --------------------------------------------------
    ch_sizes = Channel.value(file(params.chrom_sizes, checkIfExists: true))
    BIGWIG_SMOOTH(REPLISCOPE_PROFILE.out.smooth.map { comp, bed -> bed }, ch_sizes)
    BIGWIG_TREP(REPLISCOPE_TREP.out.bed, ch_sizes)

    workflow.onComplete = {
        log.info "Muller-replication-timing ${workflow.success ? 'completed' : 'FAILED'} " +
                 "in ${workflow.duration}; results in ${params.outdir}"
    }
}
