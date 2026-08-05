process ARBOR_DASHBOARD {
    tag "arbor"
    label 'process_single'

    // build_dashboard.py uses only the Python standard library (base64, gzip, json, os, sys)
    // so no container or conda env is needed — system Python 3 suffices on all HPC environments.
    // Running outside a container also avoids Singularity image pull failures on air-gapped nodes.

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
