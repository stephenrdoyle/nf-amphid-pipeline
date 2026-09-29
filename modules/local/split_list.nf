process SPLIT_LIST {
    tag "${meta.run_id}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' ?
        'docker://python:3.11-slim' : 'python:3.11-slim' }"

    input:
    tuple val(meta), path(reads)
    path assay

    output:
    // index<TAB>path, one line per input FASTQ. The index is carried so
    // a shard names its output exactly as a single task splitting the
    // whole run would have.
    tuple val(meta), path("split/read_files.txt"), emit: list
    path "split/file_report.json"                , emit: report
    path "versions.yml"                          , emit: versions

    script:
    // Which files constitute a nanopore run is a decision with rules -
    // fastq_pass in, fastq_fail out, pod5 directories ignored, MinKNOW
    // per-barcode folders walked - and they live in the tool, in
    // resolve_read_files(). Asking the tool rather than globbing in
    // Groovy keeps one answer to "what is the input". A second
    // implementation here would agree until somebody reorganised a run
    // directory, and then disagree silently.
    """
    split_concatemers.py \\
        --reads ${reads} \\
        --assay ${assay} \\
        --list-only \\
        --outdir split

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    // Real, like PREFLIGHT's stub and for the same reason. This step
    // only walks a directory and writes a list of filenames - there is
    // no work here to skip. Faking it would emit one invented path, so
    // the fan-out below it would be a fan-out of ONE in every stub
    // run, and the one thing a stub of this subworkflow is worth doing
    // is proving that N input files become N tasks and gather back.
    """
    split_concatemers.py \\
        --reads ${reads} \\
        --assay ${assay} \\
        --list-only \\
        --outdir split

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
