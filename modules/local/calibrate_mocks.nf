process CALIBRATE_MOCKS {
    tag "${run_id}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' ?
        'docker://python:3.11-slim' : 'python:3.11-slim' }"

    input:
    val  run_id
    path mocks
    path(counts, stageAs: 'assign/*')
    path reference
    // For file_resolution.tsv only. The counts tables are keyed on the
    // filename-safe sample id, because that is what reached assignment;
    // --mocks names the declared one. Without the mapping the join finds
    // nothing and the only symptom is "no mock sample has counts", which
    // reads as a bench problem rather than an id problem.
    path preflight_dir

    output:
    path "calibration/corr_factors.tsv" , emit: corr_factors
    path "calibration/mock_recovery.tsv", emit: recovery
    path "calibration/calibration.json" , emit: json
    path "reference.calibrated.tsv"     , emit: reference, optional: true
    path "versions.yml"                 , emit: versions

    script:
    """
    calibrate_mocks.py \\
        --mocks ${mocks} \\
        --assign-dir assign \\
        --reference ${reference} \\
        --preflight-dir ${preflight_dir} \\
        --run-id ${run_id} \\
        --slope-tolerance ${params.slope_tolerance} \\
        --update-reference reference.calibrated.tsv \\
        --outdir calibration

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p calibration
    touch calibration/corr_factors.tsv calibration/mock_recovery.tsv
    echo '{}' > calibration/calibration.json
    touch versions.yml
    """
}
