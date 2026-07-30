process ARBOR_DASHBOARD {
    tag "arbor"
    label 'process_single'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    path(result_files, stageAs: 'inputs/*')   // collected run outputs (flat)
    path(template)                            // dashboard/arbor_dashboard.html

    output:
    path("arbor_dashboard_loaded.html")                                                   , emit: html
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/^Python //'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // output name differs from the template (arbor_dashboard.html) to avoid a
    // staged-input vs output filename collision in the work dir
    """
    build_dashboard.py inputs ${template} arbor_dashboard_loaded.html
    """

    stub:
    """
    touch arbor_dashboard_loaded.html
    """
}
