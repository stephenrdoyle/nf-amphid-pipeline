process COMPARE_RUNS {
    tag "${label_a} vs ${label_b}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' ?
        'docker://python:3.11-slim' : 'python:3.11-slim' }"

    input:
    tuple val(label_a), path(calls_a), val(label_b), path(calls_b)

    output:
    path "comparison/comparison.json" , emit: json
    path "comparison/paired_calls.tsv", emit: paired
    path "versions.yml"               , emit: versions

    script:
    def strip = params.replicate_suffix ? "--strip '${params.replicate_suffix}'" : ''
    """
    compare_runs.py \\
        --a ${calls_a} \\
        --b ${calls_b} \\
        --label-a ${label_a} \\
        --label-b ${label_b} \\
        ${strip} \\
        --outdir comparison

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p comparison
    echo '{}' > comparison/comparison.json
    touch comparison/paired_calls.tsv
    touch versions.yml
    """
}
