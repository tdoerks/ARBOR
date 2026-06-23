process ARBOR_DASHBOARD {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.9--1'
        : 'biocontainers/python:3.9--1'}"

    input:
    path result_files, stageAs: 'results/*'
    path dashboard_template, stageAs: 'arbor_dashboard_template.html'
    path metadata

    output:
    path "arbor_dashboard.html", emit: report
    tuple val("${task.process}"), val("python"), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    // template is staged under a distinct name so the output never collides with / clobbers it
    def meta_arg = metadata ? "--metadata ${metadata}" : ''
    """
    fill_dashboard.py \\
        --template ${dashboard_template} \\
        ${meta_arg} \\
        --out arbor_dashboard.html \\
        results/
    """

    stub:
    """
    touch arbor_dashboard.html
    """
}
