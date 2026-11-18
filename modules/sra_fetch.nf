process SRA_FETCH {
    tag   "${meta.sample} (${meta.run})"
    label 'process_low'
    storeDir params.fastq_cache

    input:
    val meta

    output:
    tuple val(meta), path("${meta.sample}*.fastq.gz"), emit: reads
    path "${meta.sample}.fetch.txt",                   emit: stats

    script:
    // Reads are pulled from ENA rather than through sra-tools. ENA mirrors every SRA
    // run, serves plain FASTQ over HTTPS and is built for bulk transfer; sra-tools
    // additionally requires its SDL resolver and an S3 endpoint to be reachable.
    //
    // ENA serves the SOLiD runs in their archived colourspace form (primer base plus
    // colour calls), so those are decoded to base space on the fly - see
    // bin/decode_colourspace.py and tests/test_decode_colourspace.sh.
    def decode = meta.platform == 'ABI_SOLID'
    def pairFlag = meta.layout == 'paired' ? '--paired' : ''
    """
    ena_fetch.py --run ${meta.run} --sample ${meta.sample} ${pairFlag} \\
        --cache-dir ${params.download_cache} \\
        ${decode ? '--decode-colourspace --min-length ' + params.min_read_length : ''} \\
        --threads ${task.cpus} --stats ${meta.sample}.fetch.txt
    """
}
