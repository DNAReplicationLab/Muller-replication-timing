process REPLISCOPE_QC {
    tag   "${meta.sample}"
    label 'process_low'
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    tuple val(meta), path(counts)

    output:
    path "${meta.sample}.qc.txt", emit: stats
    path "${meta.sample}.qc.pdf", emit: plot

    script:
    """
    repliscope_qc.R ${counts} ${meta.sample}
    """
}

process REPLISCOPE_PROFILE {
    tag   "${comp.comparison}"
    label 'process_low'
    publishDir "${params.outdir}/profiles", mode: 'copy'

    input:
    tuple val(comp), path(repBed), path(nonRepBed)

    output:
    tuple val(comp), path("${comp.comparison}.ratio.tsv"),    emit: ratio
    tuple val(comp), path("${comp.comparison}.smooth.bed"),   emit: smooth
    tuple val(comp), path("${comp.comparison}.ratio.bed"),    emit: raw
    path "${comp.comparison}.normstats.txt",                  emit: stats
    path "${comp.comparison}.profile.pdf",                    emit: plot

    script:
    def value = comp.norm_value ? "--value ${comp.norm_value}" : ""
    """
    repliscope_profile.R \\
        --rep ${repBed} --nonRep ${nonRepBed} \\
        --repName ${comp.rep} --nonRepName ${comp.nonRep} \\
        --prefix ${comp.comparison} \\
        --norm ${comp.norm_mode} ${value} \\
        --loLim ${params.bin_lower_limit}
    """
}

process REPLISCOPE_TREP {
    tag   "sync-seq Trep"
    label 'process_medium'
    publishDir "${params.outdir}/profiles", mode: 'copy'

    input:
    path ratios
    val  times

    output:
    path "syncseq.trep.tsv", emit: tsv
    path "syncseq.trep.bed", emit: bed
    path "syncseq.trep.pdf", emit: plot

    script:
    """
    repliscope_trep.R \\
        --ratios ${ratios.join(',')} \\
        --times ${times.join(',')} \\
        --prefix syncseq \\
        --validate ${params.trep_validate_chrom}
    """
}
