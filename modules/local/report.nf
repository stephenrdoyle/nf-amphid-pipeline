process REPORT {
    tag "${run_id}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' ?
        'docker://python:3.11-slim' : 'python:3.11-slim' }"

    input:
    val  run_id
    path(synthesis, stageAs: 'synthesis/*')
    path comparison

    output:
    path "${run_id}_report.html", emit: report
    path "versions.yml"         , emit: versions

    script:
    def cmp = !comparison.name.startsWith('NO_') ? "--comparison ${comparison}" : ''
    """
    report.py \\
        --synthesis-dir synthesis \\
        --run-id ${run_id} \\
        ${cmp} \\
        --min-composition-reads ${params.min_composition_reads} \\
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
