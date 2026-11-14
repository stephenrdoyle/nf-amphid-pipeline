process SYNTHESISE_DB {
    tag "${run_id}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' ?
        'docker://python:3.11-slim' : 'python:3.11-slim' }"

    input:
    val  run_id
    path(clusters, stageAs: 'assign_db/*')
    path(summaries, stageAs: 'assign_db/*')
    // Contract A: the platform adapter's accounting of the run
    // (DEMUX_ACCOUNTING, docs/CONTRACTS.md), whatever the platform.
    path accounting
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
    // The protocol's hooks give the defaults of these two; passed only
    // when set.
    def host_arg = params.db_host_clade != null ? "--host-clade '${params.db_host_clade}'" : ''
    def rel_arg  = params.db_min_rel_abundance != null ? "--min-rel-abundance ${params.db_min_rel_abundance}" : ''
    // The demultiplexer's index matrix: unused index pairs measure
    // misassignment directly, a cross-check on the fitted rates.
    def matrix_arg = !index_matrix.name.startsWith('NO_') ? "--index-matrix ${index_matrix}" : ''
    """
    synthesise_db.py \\
        --assign-dir assign_db \\
        --preflight-dir ${preflight_dir} \\
        --database ${database} \\
        --accounting ${accounting} ${matrix_arg} \\
        --protocol ${params.db_protocol} \\
        --calling-rule ${params.db_calling_rule} \\
        --run-id ${run_id} \\
        --min-depth ${params.min_depth} \\
        --absolute-floor ${params.absolute_floor} \\
        ${rel_arg} ${host_arg} \\
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
