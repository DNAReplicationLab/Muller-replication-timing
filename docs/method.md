# Method

## 1. Acquisition (`modules/sra_fetch.nf`)

Runs are fetched with `prefetch` and converted to FASTQ. Illumina runs use `fasterq-dump`;
SOLiD runs use `fastq-dump -B`, which decodes colourspace to base space (see
[deviations](deviations.md)). Output is cached in `params.fastq_cache` with Nextflow's
`storeDir`, so re-running the workflow does not re-download.

## 2. Alignment (`modules/bowtie2_align.nf`)

Bowtie 2 against the reference index, piped straight into `samtools sort`. Paired-end data
additionally go through `samtools fixmate -m` and `samtools markdup`, so PCR and optical
duplicates are marked. Single-end data are not deduplicated: for single-end reads a
duplicate cannot be distinguished from an independent fragment with the same 5′ end, and at
these depths that is common by chance.

## 3. Counting (`modules/count_bins.nf`)

This is the step that defines the measurement. For each sample:

1. `bedtools makewindows` divides the genome into fixed bins (1 kb by default).
2. Alignments are filtered to MAPQ ≥ 30, `-F 3840` (no secondary, QC-fail, duplicate or
   supplementary alignments), and for paired data `-f 64` (first in pair only, so each
   fragment is counted once). Any read carrying an `XS:i:` tag — meaning Bowtie 2 found a
   second alignment for it — is removed, leaving only uniquely placed reads.
3. `bedtools genomecov -5` reduces each surviving read to its 5′ end position.
4. `bedtools map ... -o sum` sums those positions into the bins.

The result is a BED5 file: chromosome, bin start, bin end, sample name, read count.
This is exactly the input `Repliscope::loadBed` expects.

Counting 5′ ends rather than coverage matters: coverage would weight a bin by read length
and by how reads straddle bin boundaries, whereas a 5′ end is one unambiguous observation
of one fragment.

## 4. Profiles (`bin/repliscope_profile.R`)

Per the protocol, step 61A:

1. `loadBed` with the sacCer3 mask, restricting to the nuclear chromosomes.
2. `rmOutliers(method="median", loLim=0.25)` drops bins holding less than a quarter of the
   genome-wide median count — repeats, low-mappability regions and the rDNA, which give
   noisy ratios. This is the same criterion as the 2014 paper's "windows where fewer than a
   quarter of the expected number of reads were mapped in either sample were excluded".
3. `makeRatio` divides replicating by non-replicating counts per bin and corrects for the
   difference in total reads.
4. `normaliseRatio` puts the ratio onto the relative copy number scale — see the README for
   the three modes.
5. `smoothRatio` fits a cubic spline within each run of adjacent populated bins.

## 5. Median replication time (`bin/repliscope_trep.R`)

For a sync-seq time course, each bin's relative copy number is followed across timepoints
and a sigmoid is fitted, anchored at copy number 1 at t = 0. Trep is the time at which that
bin reached copy number 1.5 — replicated in half the cells. Bins whose fitted Trep falls
outside the sampled window are discarded as unsupported.

## 6. Tracks (`modules/bigwig.nf`)

Smoothed profiles and the Trep track are written as bigWig against the reference
`chrom.sizes`, ready for a UCSC track hub or JBrowse.
