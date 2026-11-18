process BOWTIE2_ALIGN {
    tag   "${meta.sample}"
    label 'process_high'
    // 'link' not 'copy': the BAMs run to tens of GB and a hard link costs nothing
    publishDir "${params.outdir}/bam",  mode: 'link', pattern: '*.{bam,bai}'
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.{log,txt}'

    input:
    tuple val(meta), path(reads)
    tuple path(index), val(index_base)

    output:
    tuple val(meta), path("${meta.sample}.bam"), path("${meta.sample}.bam.bai"), emit: bam
    path  "${meta.sample}.bowtie2.log",                                          emit: log
    path  "${meta.sample}.markdup.txt",                  optional: true,         emit: metrics

    script:
    // localMapper.sh uses `bowtie2 -a` for single-end reads (report every alignment found).
    // Omitted here: multi-mapping reads are discarded downstream by the MAPQ>=30 and XS:i:
    // filters, so -a changes nothing in the counts while costing orders of magnitude in
    // runtime and disk. Duplicate marking uses samtools markdup rather than Picard
    // MarkDuplicates. Both are recorded in docs/deviations.md.
    if (meta.layout == 'paired')
        """
        bowtie2 -p ${task.cpus} -q --phred33 -X ${params.max_insert} \\
            -x ${index_base} -1 ${reads[0]} -2 ${reads[1]} 2> ${meta.sample}.bowtie2.log \\
          | samtools sort -n -@ ${task.cpus} -m ${params.sort_memory} -T nsort -O bam - \\
          | samtools fixmate -m -@ ${task.cpus} - - \\
          | samtools sort -@ ${task.cpus} -m ${params.sort_memory} -T csort -O bam - \\
          | samtools markdup -@ ${task.cpus} -f ${meta.sample}.markdup.txt - ${meta.sample}.bam
        samtools index -@ ${task.cpus} ${meta.sample}.bam
        """
    else
        """
        bowtie2 -p ${task.cpus} -q --phred33 \\
            -x ${index_base} -U ${reads} 2> ${meta.sample}.bowtie2.log \\
          | samtools sort -@ ${task.cpus} -m ${params.sort_memory} -T csort -O bam -o ${meta.sample}.bam -
        samtools index -@ ${task.cpus} ${meta.sample}.bam
        """
}
