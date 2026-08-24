process CAT_FASTQ {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(r1_files), path(r2_files)

    output:
    tuple val(meta), path("${meta.id}_R1.fastq.gz"), path("${meta.id}_R2.fastq.gz"), emit: reads

    script:
    // gzip files can be concatenated directly — result is a valid gzip stream
    """
    cat ${r1_files.join(' ')} > ${meta.id}_R1.fastq.gz
    cat ${r2_files.join(' ')} > ${meta.id}_R2.fastq.gz
    """
}
