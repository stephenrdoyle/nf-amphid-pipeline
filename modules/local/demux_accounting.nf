process DEMUX_ACCOUNTING {
    tag "${run_id}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' ?
        'docker://python:3.11-slim' : 'python:3.11-slim' }"

    input:
    val  run_id
    path demux_summary
    path split_summary
    path preflight_dir

    output:
    path "demux/demux_accounting.json", emit: json
    path "versions.yml"               , emit: versions

    script:
    // Contract A's accounting (docs/CONTRACTS.md): where the run's reads
    // went before any sample saw them, in fates every platform shares,
    // and each well's plate position. Nanopore's demux and split
    // summaries are folded in; an Illumina or predemultiplexed run has
    // neither (their NO_ placeholders), and its accounting says that
    // nothing is known of reads that reached no sample.
    def demux_arg = !demux_summary.name.startsWith('NO_') ? "--demux-summary ${demux_summary}" : ''
    def split_arg = !split_summary.name.startsWith('NO_') ? "--split-summary ${split_summary}" : ''
    """
    mkdir -p demux
    demux_accounting.py \\
        --run-id ${run_id} \\
        --preflight-dir ${preflight_dir} \\
        ${demux_arg} ${split_arg} \\
        --out demux/demux_accounting.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p demux
    echo '{"run_id":"${run_id}"}' > demux/demux_accounting.json
    touch versions.yml
    """
}
