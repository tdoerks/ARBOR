process SPADES_ASSEMBLE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/spades:4.0.0--h5ca1c30_1' :
        'quay.io/biocontainers/spades:4.0.0--h5ca1c30_1' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}.scaffolds.fasta"), optional: true, emit: scaffolds
    tuple val(meta), path("${prefix}.assembly_stats.tsv"),              emit: stats
    tuple val("${task.process}"), val('spades'), eval("spades.py --version 2>&1 | sed 's/SPAdes v//'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    prefix    = task.ext.prefix ?: "${meta.id}"
    def reads_args = meta.single_end
        ? "-s ${reads}"
        : "-1 ${reads[0]} -2 ${reads[1]}"
    def mem_gb = (task.memory.toGiga() as int)
    """
    spades.py \\
        --careful \\
        $reads_args \\
        -o spades_out \\
        -t $task.cpus \\
        -m $mem_gb \\
        $args

    if [ -f spades_out/scaffolds.fasta ]; then
        cp spades_out/scaffolds.fasta ${prefix}.scaffolds.fasta
    else
        touch ${prefix}.scaffolds.fasta
    fi

    awk -v sample="${meta.id}" '
    BEGIN { n=0; total=0; largest=0 }
    /^>/ {
        if (seq != "") {
            len = length(seq); n++; total += len
            if (len > largest) largest = len
            lens[n] = len
        }
        seq = ""; next
    }
    { seq = seq \$0 }
    END {
        if (seq != "") {
            len = length(seq); n++; total += len
            if (len > largest) largest = len; lens[n] = len
        }
        n50 = 0
        if (n > 0) {
            asort(lens)
            cumsum = 0
            for (i = n; i >= 1; i--) {
                cumsum += lens[i]
                if (cumsum >= total / 2) { n50 = lens[i]; break }
            }
        }
        print "sample\\tn_contigs\\ttotal_length\\tlargest_contig\\tN50"
        print sample "\\t" n "\\t" total "\\t" largest "\\t" n50
    }' ${prefix}.scaffolds.fasta > ${prefix}.assembly_stats.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.scaffolds.fasta
    printf 'sample\\tn_contigs\\ttotal_length\\tlargest_contig\\tN50\\n${meta.id}\\t0\\t0\\t0\\t0\\n' > ${prefix}.assembly_stats.tsv
    """
}
