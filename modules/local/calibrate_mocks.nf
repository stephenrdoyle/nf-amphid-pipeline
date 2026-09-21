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
