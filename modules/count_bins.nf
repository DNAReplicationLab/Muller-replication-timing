process COUNT_BINS {
    tag   "${meta.sample} ${params.bin_size}bp"
    label 'process_medium'
    publishDir "${params.outdir}/counts", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.sample}.${params.bin_size}bp.bed"), emit: counts
    path "${meta.sample}.${params.bin_size}bp.stats.txt",             emit: stats

    script:
    // Reproduces the counting logic of localMapper.sh:
    //   -q 30         unique, confidently placed alignments
    //   -F 3840       drop secondary, QC-fail, duplicate and supplementary alignments
    //   -f 64         first in pair only (paired data), so each fragment counts once
    //   grep -v XS:i: drop any read with a reported secondary alignment
    //   genomecov -5  count the 5' end position of each surviving read
    def pairflag = meta.layout == 'paired' ? '-f 64' : ''
    """
    samtools idxstats ${bam} | awk 'BEGIN{OFS="\\t"} \$2>0 {print \$1,\$2}' > genome.txt
    bedtools makewindows -g genome.txt -w ${params.bin_size} > windows.bed

    samtools view -h -@ ${task.cpus} -q 30 -F 3840 ${pairflag} -L windows.bed ${bam} \\
        | grep -v XS:i: \\
        | samtools view -@ ${task.cpus} -b - \\
        | bedtools genomecov -5 -d -ibam stdin \\
        | awk 'BEGIN{OFS="\\t"} \$3>0 {print \$1,\$2-1,\$2,"5prime",\$3}' > fiveprime.bed

    bedtools map -a windows.bed -b fiveprime.bed -c 5 -o sum -null 0 \\
        | awk 'BEGIN{OFS="\\t"} \$4>0 {print \$1,\$2,\$3,"${meta.sample}",\$4}' \\
        > ${meta.sample}.${params.bin_size}bp.bed

    {
      echo -e "sample\\t${meta.sample}"
      echo -e "run\\t${meta.run}"
      echo -e "platform\\t${meta.platform}"
      echo -e "total_alignments\\t\$(samtools view -c ${bam})"
      echo -e "counted_5prime_ends\\t\$(awk '{s+=\$5} END {print s+0}' fiveprime.bed)"
      echo -e "binned_reads\\t\$(awk '{s+=\$5} END {print s+0}' ${meta.sample}.${params.bin_size}bp.bed)"
      echo -e "nonzero_bins\\t\$(wc -l < ${meta.sample}.${params.bin_size}bp.bed)"
      echo -e "total_bins\\t\$(wc -l < windows.bed)"
    } > ${meta.sample}.${params.bin_size}bp.stats.txt
    """
}
