process REPORT_DB {
    tag "${run_id}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' ?
        'docker://python:3.11-slim' : 'python:3.11-slim' }"

    input:
    val  run_id
    path(synthesis, stageAs: 'synthesis/*')

    output:
    path "${run_id}_report.html", emit: report
    path "versions.yml"         , emit: versions

    script:
    """
    report_db.py \\
        --synthesis-dir synthesis \\
        --run-id ${run_id} \\
        --out ${run_id}_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    echo '<html></html>' > ${run_id}_report.html
    touch versions.yml
    """
}
