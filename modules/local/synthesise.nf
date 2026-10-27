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
    path split_summary
    path samplesheet
    path reference
    path preflight_dir

    output:
    path "synthesis/calls.tsv"      , emit: calls
    path "synthesis/calibration.tsv", emit: calibration
    path "synthesis/read_fates.tsv" , emit: fates
    // Optional because it is only written when there is a grouping to
    // write - a metadata file, or a subject with more than one library.
    // An empty file would say "no groups" and "no metadata" with the
    // same bytes.
    path "synthesis/group_summary.tsv", emit: groups, optional: true
    path "synthesis/run_qc.json"    , emit: qc
    path "versions.yml"             , emit: versions

    script:
    // The funnel starts wherever the platform's first count is: preparation
    // on Illumina, demultiplexing on nanopore. Passing only what exists
    // means the boundary check is made where it can be, and reported as
    // not made where it cannot.
    def prep_arg  = preps ? "--prep-dir prep" : ''
    def demux_arg = !demux_summary.name.startsWith('NO_') ? "--demux-summary ${demux_summary}" : ''
    // Without this the funnel starts at post-split amplicon units rather
    // than at reads off the sequencer, and every percentage below it is
    // taken against a denominator the splitting quietly changed.
    def split_arg = !split_summary.name.startsWith('NO_') ? "--split-summary ${split_summary}" : ''
    // The resolved metadata join already sits inside the staged
    // pre-flight directory, so it needs no channel of its own. Passed
    // by path rather than derived inside synthesise.py, because a tool
    // that goes looking inside a directory for optional inputs is a
    // tool whose inputs are not declared anywhere.
    //
    // --slope-tolerance is params.slope_tolerance, documented as the
    // tolerated deviation in the standard-curve slope - which is this
    // process's slope. Without it synthesise.py fell back on a default of
    // its own and the parameter silently did nothing here.
    """
    synthesise.py \\
        --assign-dir assign \\
        --samplesheet ${samplesheet} \\
        --reference ${reference} \\
        --preflight-dir ${preflight_dir} \\
        --sample-metadata ${preflight_dir}/sample_metadata.tsv \\
        ${prep_arg} ${demux_arg} ${split_arg} \\
        --run-id ${run_id} \\
        ${params.spiked.toString().toLowerCase() == 'true' ? '--spiked' : ''} \\
        --min-depth ${params.min_depth} \\
        --absolute-floor ${params.absolute_floor} \\
        --slope-tolerance ${params.slope_tolerance} \\
        --permutations ${params.permutations} \\
        --plasmid-warn ${params.plasmid_warn} \\
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
    touch synthesis/group_summary.tsv
    echo '{"run_id":"${run_id}"}' > synthesis/run_qc.json
    touch versions.yml
    """
}
