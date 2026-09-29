process MERGE_SPLIT {
    tag "${meta.run_id}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' ?
        'docker://python:3.11-slim' : 'python:3.11-slim' }"

    input:
    val meta
    path parts, stageAs: 'parts/*'
    path file_report
    val n_expected

    output:
    path "split/split_summary.json", emit: summary
    path "versions.yml"            , emit: versions

    script:
    // --expect-parts is the point of this step as much as the summing
    // is. Every percentage downstream is taken against reads_in, so a
    // shard that quietly failed to arrive would not raise an error -
    // it would produce a read funnel that balances perfectly against a
    // smaller run than the one that was sequenced.
    """
    merge_split.py \\
        --parts parts \\
        --file-report ${file_report} \\
        --min-length ${params.min_unit_length} \\
        --anchor-max-edits ${params.anchor_max_edits} \\
        --expect-parts ${n_expected} \\
        --outdir split

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p split
    echo '{"reads_in":0,"pieces_out":0,"parts":${n_expected}}' > split/split_summary.json
    touch versions.yml
    """
}
