process COMPARE_RUNS {
    tag "${label_a} vs ${label_b}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' ?
        'docker://python:3.11-slim' : 'python:3.11-slim' }"

    input:
    // Each side staged into its own directory: both are calls.tsv, and
    // two inputs of one name collide when Nextflow stages them.
    tuple val(label_a), path(calls_a, stageAs: 'a/*'), val(label_b), path(calls_b, stageAs: 'b/*')

    output:
    path "comparison/comparison.json" , emit: json
    path "comparison/paired_calls.tsv", emit: paired
    path "versions.yml"               , emit: versions

    script:
    def strip = params.replicate_suffix ? "--strip '${params.replicate_suffix}'" : ''
    // An empty label is left off, so the tool takes it from the run_id
    // column of that side's calls.tsv.
    def lab_a = label_a ? "--label-a '${label_a}'" : ''
    def lab_b = label_b ? "--label-b '${label_b}'" : ''
    """
    compare_runs.py \\
        --a ${calls_a} \\
        --b ${calls_b} \\
        ${lab_a} ${lab_b} \\
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
