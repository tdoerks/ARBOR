process SPADES_ASSEMBLE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/spades:4.0.0--haf24da9_4' :
        'quay.io/biocontainers/spades:4.0.0--haf24da9_4' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}.scaffolds.fasta"), optional: true, emit: scaffolds
    tuple val(meta), path("${prefix}.assembly_stats.tsv"),              emit: stats
    tuple val("${task.process}"), val('spades'), eval("spades.py --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+'"), topic: versions

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

    python3 - "${meta.id}" ${prefix}.scaffolds.fasta ${prefix}.assembly_stats.tsv << 'PYEOF'
import sys
T = chr(9); NL = chr(10)
sample, fasta, out = sys.argv[1], sys.argv[2], sys.argv[3]
lens = []
seq = []
with open(fasta) as fh:
    for line in fh:
        line = line.rstrip()
        if line.startswith('>'):
            if seq:
                lens.append(len(''.join(seq)))
            seq = []
        else:
            seq.append(line)
    if seq:
        lens.append(len(''.join(seq)))
n = len(lens)
total = sum(lens)
largest = max(lens) if lens else 0
n50 = 0
if lens:
    lens_sorted = sorted(lens, reverse=True)
    cumsum = 0
    for l in lens_sorted:
        cumsum += l
        if cumsum >= total / 2:
            n50 = l; break
with open(out, 'w') as fh:
    fh.write(T.join(['sample','n_contigs','total_length','largest_contig','N50']) + NL)
    fh.write(T.join([sample, str(n), str(total), str(largest), str(n50)]) + NL)
PYEOF
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.scaffolds.fasta
    printf 'sample\\tn_contigs\\ttotal_length\\tlargest_contig\\tN50\\n${meta.id}\\t0\\t0\\t0\\t0\\n' > ${prefix}.assembly_stats.tsv
    """
}
