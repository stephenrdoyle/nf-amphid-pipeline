process PREFLIGHT {
    tag "preflight"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' ?
        'docker://python:3.11-slim' : 'python:3.11-slim' }"

    input:
    path reference
    path samplesheet
    path runs
    path index_map
    path assay
    path primers
    path mocks

    output:
    path "preflight"                          , emit: dir
    path "preflight/reference.fasta"          , emit: fasta
    // The primer-trimmed panel: reads arrive with their primers removed,
    // so this is what assignment competes them against. Optional because
    // it needs --primers, and a panel can be validated without them.
    path "preflight/reference.insert.fasta"   , emit: insert_fasta, optional: true
    path "preflight/tag_report.tsv"           , emit: tag_report
    path "preflight/tag_report.txt"           , emit: tag_report_txt
    path "preflight/file_resolution.tsv"      , emit: file_resolution
    path "preflight/preflight_summary.json"   , emit: summary
    path "preflight/validation.log"           , emit: log
    path "versions.yml"                       , emit: versions

    script:
    // The run manifest is optional at the interface level: --input_dir builds
    // a single-row one. Nextflow needs a placeholder file either way, so an
    // absent manifest arrives as the empty file assets/NO_FILE.
    def manifest_arg = runs.name != 'NO_FILE'
        ? "--runs ${runs}"
        : "--input-dir ${params.input_dir} --run-id ${params.run_id} " +
          "--platform ${params.platform}"
    // Overrides are omitted entirely when unset, so the assay config decides.
    // An explicitly empty expected_tag means "skip the cross-check", which is
    // different from "not specified" - hence the null test rather than truthiness.
    def tag_edits = params.tag_max_edits != null ? "--tag-max-edits ${params.tag_max_edits}" : ''
    def exp_tag   = params.expected_tag  != null ? "--expected-tag '${params.expected_tag}'" : ''
    def reps      = params.barcode_repeats != null ? "--barcode-repeats ${params.barcode_repeats}" : ''
    def prim      = primers.name != 'NO_FILE' ? "--primers ${primers}" : ''
    def mock_arg  = mocks.name   != 'NO_FILE' ? "--mocks ${mocks}" : ''
    def read_len  = params.illumina_read_length != null
        ? "--illumina-read-length ${params.illumina_read_length}" : ''
    """
    preflight.py \\
        --reference ${reference} \\
        --samplesheet ${samplesheet} \\
        ${manifest_arg} \\
        --index-map ${index_map} \\
        --assay ${assay} \\
        ${prim} ${mock_arg} ${read_len} \\
        ${tag_edits} ${exp_tag} ${reps} \\
        --outdir preflight

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p preflight
    touch preflight/reference.fasta preflight/tag_report.tsv preflight/tag_report.txt
    touch preflight/file_resolution.tsv preflight/validation.log
    touch preflight/reference.insert.fasta
    echo '{}' > preflight/preflight_summary.json
    touch versions.yml
    """
}
