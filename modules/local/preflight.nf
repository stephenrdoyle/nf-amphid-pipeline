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
    path metadata

    output:
    path "preflight"                          , emit: dir
    path "preflight/reference.fasta"          , emit: fasta
    // The primer-trimmed panel: reads arrive with their primers removed,
    // so this is what assignment competes them against. NOT optional.
    // It used to be, and a missing one left ASSIGN with an empty input,
    // so every sample skipped assignment and the run printed Complete.
    // main.nf now requires --primers, and preflight.py fails when the
    // primers locate nowhere, so its absence here is a bug worth a stop.
    path "preflight/reference.insert.fasta"   , emit: insert_fasta
    path "preflight/tag_report.tsv"           , emit: tag_report
    path "preflight/tag_report.txt"           , emit: tag_report_txt
    path "preflight/runs_resolved.tsv" , emit: runs_resolved
    path "preflight/file_resolution.tsv"      , emit: file_resolution
    // The grouping join, resolved to one row per library. Written even
    // with no --metadata, because subject_id is a grouping field the
    // samplesheet already declares. Optional only because a
    // reference-only run has no samplesheet to join to.
    path "preflight/sample_metadata.tsv"      , emit: sample_metadata, optional: true
    path "preflight/preflight_summary.json"   , emit: summary
    path "preflight/validation.log"           , emit: log
    path "versions.yml"                       , emit: versions

    script:
    // The run manifest is optional at the interface level: --input_dir builds
    // a single-row one. Nextflow needs a placeholder file either way, so an
    // absent manifest arrives as the empty file assets/NO_FILE.
    // Nextflow stages the manifest, so by the time pre-flight reads it
    // the file is in a scratch directory and any relative reads path in
    // it would resolve there. The original directory is known here and
    // nowhere else, so it is passed explicitly.
    // A relative --input_dir means relative to where the pipeline was
    // launched, but this script runs in a work/ scratch directory, so the
    // launch directory is passed as --path-base and preflight.py makes the
    // path absolute against it before writing it into runs.derived.tsv.
    // (Without it the path ended up resolved against <work>/preflight/.)
    // The resolution itself stays in Python, where it is tested; Groovy
    // only says where "here" was. Paths are single-quoted because a
    // directory name with a space would otherwise split into two arguments.
    def base_arg = params.runs
        ? "--path-base '${file(params.runs).parent}'"
        : (params.input_dir ? "--path-base '${workflow.launchDir}'" : '')
    def manifest_arg = !runs.name.startsWith('NO_')
        ? "--runs ${runs}"
        : "--input-dir '${params.input_dir}' --run-id '${params.run_id}' " +
          "--platform ${params.platform}"
    // Overrides are omitted entirely when unset, so the assay config decides.
    // An explicitly empty expected_tag means "skip the cross-check", which is
    // different from "not specified" - hence the null test rather than truthiness.
    def tag_edits = params.tag_max_edits != null ? "--tag-max-edits ${params.tag_max_edits}" : ''
    def exp_tag   = params.expected_tag  != null ? "--expected-tag '${params.expected_tag}'" : ''
    def reps      = params.barcode_repeats != null ? "--barcode-repeats ${params.barcode_repeats}" : ''
    def prim      = !primers.name.startsWith('NO_') ? "--primers ${primers}" : ''
    def mock_arg  = !mocks.name.startsWith('NO_') ? "--mocks ${mocks}" : ''
    def meta_arg  = !metadata.name.startsWith('NO_') ? "--metadata ${metadata}" : ''
    def read_len  = params.illumina_read_length != null
        ? "--illumina-read-length ${params.illumina_read_length}" : ''
    def idx_edit  = params.min_index_edit != null
        ? "--min-index-edit ${params.min_index_edit}" : ''
    """
    preflight.py \\
        --reference ${reference} \\
        --samplesheet ${samplesheet} \\
        ${manifest_arg} \\
        --index-map ${index_map} \\
        --assay ${assay} \\
        ${params.spiked.toString().toLowerCase() == 'true' ? '--spiked' : ''} \\
        ${prim} ${mock_arg} ${meta_arg} ${base_arg} ${read_len} \\
        ${idx_edit} \\
        ${tag_edits} ${exp_tag} ${reps} \\
        --outdir preflight

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    // The stub RUNS PRE-FLIGHT FOR REAL, and that is deliberate.
    //
    // A stub does not read its inputs, so a hardcoded file_resolution
    // emits whatever rows were written into this file - which meant
    // `-profile test_ont -stub-run` produced an Illumina sample called
    // STUB_SAMPLE and exercised the Illumina branch, whichever profile
    // was asked for. A stub that cannot reflect its own inputs cannot
    // test a workflow whose routing IS its inputs: the platform column
    // decides which branch runs, so a stub that invents the platform
    // column decides it too.
    //
    // There is also nothing here worth stubbing. `-stub-run` exists to
    // skip expensive work; pre-flight is standard-library Python that
    // touches no reads and finishes in about a second. It is the one
    // process in this pipeline whose real execution is cheaper than a
    // convincing fake of it.
    // Same --path-base / quoting rule as the script block above.
    def base_arg = params.runs
        ? "--path-base '${file(params.runs).parent}'"
        : (params.input_dir ? "--path-base '${workflow.launchDir}'" : '')
    def manifest_arg = !runs.name.startsWith('NO_')
        ? "--runs ${runs}"
        : "--input-dir '${params.input_dir}' --run-id '${params.run_id}' " +
          "--platform ${params.platform}"
    def tag_edits = params.tag_max_edits != null ? "--tag-max-edits ${params.tag_max_edits}" : ''
    def exp_tag   = params.expected_tag  != null ? "--expected-tag '${params.expected_tag}'" : ''
    def reps      = params.barcode_repeats != null ? "--barcode-repeats ${params.barcode_repeats}" : ''
    def prim      = !primers.name.startsWith('NO_') ? "--primers ${primers}" : ''
    def mock_arg  = !mocks.name.startsWith('NO_') ? "--mocks ${mocks}" : ''
    def meta_arg  = !metadata.name.startsWith('NO_') ? "--metadata ${metadata}" : ''
    def read_len  = params.illumina_read_length != null
        ? "--illumina-read-length ${params.illumina_read_length}" : ''
    def idx_edit  = params.min_index_edit != null
        ? "--min-index-edit ${params.min_index_edit}" : ''
    """
    preflight.py \\
        --reference ${reference} \\
        --samplesheet ${samplesheet} \\
        ${manifest_arg} \\
        --index-map ${index_map} \\
        --assay ${assay} \\
        ${params.spiked.toString().toLowerCase() == 'true' ? '--spiked' : ''} \\
        ${prim} ${mock_arg} ${meta_arg} ${base_arg} ${read_len} \\
        ${idx_edit} \\
        ${tag_edits} ${exp_tag} ${reps} \\
        --outdir preflight

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
