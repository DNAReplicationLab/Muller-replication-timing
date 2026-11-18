process BED_TO_BIGWIG {
    tag   "${bed.simpleName}"
    label 'process_low'
    publishDir "${params.outdir}/bigwig", mode: 'copy'

    input:
    path bed
    path chrom_sizes

    output:
    path "*.bw", emit: bigwig

    script:
    """
    bed_to_bigwig.py --bed ${bed} --chrom-sizes ${chrom_sizes} \\
        --out ${bed.getName().replaceAll(/\.bed$/, '')}.bw
    """
}
