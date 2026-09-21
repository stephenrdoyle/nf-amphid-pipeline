process SYNTHESISE {
    tag "${run_id}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' ?
        'docker://python:3.11-slim' : 'python:3.11-slim' }"

    input:
    val  run_id
    path(counts, stageAs: 'assign/*')
    path(summaries, stageAs: 'assign/*')
    path(preps, stageAs: 'prep/*')
    path demux_summary
    path samplesheet
    path reference
    path preflight_dir

    output:
    path "synthesis/calls.tsv"      , emit: calls
    path "synthesis/calibration.tsv", emit: calibration
    path "synthesis/read_fates.tsv" , emit: fates
    path "synthesis/run_qc.json"    , emit: qc
    path "versions.yml"             , emit: versions

    script:
    // The funnel starts wherever the platform's first count is: preparation
    // on Illumina, demultiplexing on nanopore. Passing only what exists
    // means the boundary check is made where it can be, and reported as
    // not made where it cannot.
    def prep_arg  = preps ? "--prep-dir prep" : ''
    def demux_arg = demux_summary.name != 'NO_FILE' ? "--demux-summary ${demux_summary}" : ''
    """
    synthesise.py \\
        --assign-dir assign \\
        --samplesheet ${samplesheet} \\
        --reference ${reference} \\
        --preflight-dir ${preflight_dir} \\
        ${prep_arg} ${demux_arg} \\
        --run-id ${run_id} \\
        --min-depth ${params.min_depth} \\
        --absolute-floor ${params.absolute_floor} \\
        --permutations ${params.permutations} \\
        --outdir synthesis

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p synthesis
    touch synthesis/calls.tsv synthesis/calibration.tsv synthesis/read_fates.tsv
    echo '{"run_id":"${run_id}"}' > synthesis/run_qc.json
    touch versions.yml
    """
}
