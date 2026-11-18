process BOWTIE2_INDEX {
    tag   "${fasta.simpleName}"
    label 'process_medium'
    storeDir "${params.reference_dir}/bowtie2"

    input:
    path fasta

    output:
    tuple path("${fasta.simpleName}*.bt2"), val("${fasta.simpleName}"), emit: index

    script:
    """
    bowtie2-build --threads ${task.cpus} -f ${fasta} ${fasta.simpleName}
    """
}
