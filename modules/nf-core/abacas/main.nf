process ABACAS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/abacas:1.3.1--pl526_0' :
        'quay.io/biocontainers/abacas:1.3.1--pl526_0' }"

    input:
    tuple val(meta), path(scaffold)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("${prefix}.*"), emit: results
    tuple val("${task.process}"), val('abacas'), eval("abacas.pl --version 2>&1 | grep 'ABACAS\\.' | sed 's/ABACAS\\.//' || true"), topic: versions, emit: versions_abacas
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}.abacas"
    """
    # ABACAS requires a single-record reference — split multi-segment FASTA into one file per record
    python3 - << 'PYEOF'
segs = {}
current = None
with open("${fasta}") as f:
    for line in f:
        line = line.rstrip()
        if line.startswith(">"):
            current = line[1:].split()[0]
            segs[current] = []
        elif current:
            segs[current].append(line)
for seg_id, lines in segs.items():
    with open(f"{seg_id}.split.fa", "w") as out:
        out.write(f">{seg_id}\\n" + "\\n".join(lines) + "\\n")
PYEOF

    # Run ABACAS once per segment
    for ref_seg in *.split.fa; do
        seg=\$(basename "\$ref_seg" .split.fa)
        abacas.pl \\
            -r "\$ref_seg" \\
            -q ${scaffold} \\
            ${args} \\
            -o ${prefix}.\${seg}

        [ -f "${prefix}.\${seg}.bin" ]  && sort "${prefix}.\${seg}.bin" > "${prefix}.\${seg}.bin.tmp" \\
                                        && mv "${prefix}.\${seg}.bin.tmp" "${prefix}.\${seg}.bin" || true
        [ -f nucmer.delta ]             && mv nucmer.delta             ${prefix}.\${seg}.nucmer.delta          || true
        [ -f nucmer.filtered.delta ]    && mv nucmer.filtered.delta    ${prefix}.\${seg}.nucmer.filtered.delta || true
        [ -f nucmer.tiling ]            && mv nucmer.tiling            ${prefix}.\${seg}.nucmer.tiling         || true
        [ -f unused_contigs.out ]       && mv unused_contigs.out       ${prefix}.\${seg}.unused.contigs.out    || true
    done

    rm -f *.split.fa
    """

    stub:
    prefix   = task.ext.prefix ?: "${meta.id}.abacas"
    """
    for seg in S M L; do
        touch ${prefix}.\${seg}.fasta
        touch ${prefix}.\${seg}.crunch
        touch ${prefix}.\${seg}.gaps
        touch ${prefix}.\${seg}.gaps.tab
        touch ${prefix}.\${seg}.tab
        touch ${prefix}.\${seg}.bin
        touch ${prefix}.\${seg}.nucmer.delta
        touch ${prefix}.\${seg}.nucmer.filtered.delta
        touch ${prefix}.\${seg}.nucmer.tiling
        touch ${prefix}.\${seg}.unused.contigs.out
    done
    """
}
