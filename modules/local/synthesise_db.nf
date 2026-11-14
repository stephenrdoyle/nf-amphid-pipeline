process SYNTHESISE_DB {
    tag "${run_id}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' ?
        'docker://python:3.11-slim' : 'python:3.11-slim' }"

    input:
    val  run_id
    path(clusters, stageAs: 'assign_db/*')
    path(summaries, stageAs: 'assign_db/*')
    path demux_summary
    path split_summary
    path index_matrix
    path preflight_dir
    path database

    output:
    path "synthesis/db_calls.tsv"   , emit: calls
    path "synthesis/db_wells.tsv"   , emit: wells
    path "synthesis/db_prevalence.tsv", emit: prevalence
    path "synthesis/db_run_qc.json" , emit: qc
    path "versions.yml"             , emit: versions

    script:
    // Database mode's synthesis: relative 18S abundance at the rank each
    // call reached, judged against the run's blanks and a relative
    // abundance floor, with off-target fractions per well. No standards,
    // no plasmids - see bin/synthesise_db.py.
    def demux_arg = !demux_summary.name.startsWith('NO_') ? "--demux-summary ${demux_summary}" : ''
    def split_arg = !split_summary.name.startsWith('NO_') ? "--split-summary ${split_summary}" : ''
    // The demultiplexer's index matrix: unused index pairs measure
    // misassignment directly, a cross-check on the fitted rates.
    def matrix_arg = !index_matrix.name.startsWith('NO_') ? "--index-matrix ${index_matrix}" : ''
    """
    synthesise_db.py \\
        --assign-dir assign_db \\
        --preflight-dir ${preflight_dir} \\
        --database ${database} \\
        ${demux_arg} ${split_arg} ${matrix_arg} \\
        --run-id ${run_id} \\
        --min-depth ${params.min_depth} \\
        --absolute-floor ${params.absolute_floor} \\
        --min-rel-abundance ${params.db_min_rel_abundance} \\
        --host-clade '${params.db_host_clade}' \\
        --misassign-default-rate ${params.db_misassign_default_rate} \\
        --misassign-margin ${params.db_misassign_margin} \\
        --outdir synthesis

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p synthesis
    touch synthesis/db_calls.tsv synthesis/db_wells.tsv synthesis/db_prevalence.tsv
    echo '{"run_id":"${run_id}"}' > synthesis/db_run_qc.json
    touch versions.yml
    """
}
