#!/usr/bin/env nextflow

/*
 * nf-helminth-amplicon
 *
 * Targeted amplicon panel analysis on Illumina and nanopore, built around
 * competitive mapping against a defined reference panel with spiked
 * quantification standards and tagged control plasmids.
 *
 * Increment 1: pre-flight validation only. Every input is checked and
 * everything derivable is derived before a single read is touched.
 */

nextflow.enable.dsl = 2

include { PREFLIGHT      } from './modules/local/preflight'
include { ILLUMINA_PATH  } from './subworkflows/local/illumina'
include { NANOPORE_PATH  } from './subworkflows/local/nanopore'
include { ASSIGN         } from './modules/local/assign'
include { SYNTHESISE     } from './modules/local/synthesise'
include { CALIBRATE_MOCKS} from './modules/local/calibrate_mocks'
include { REPORT         } from './modules/local/report'

def helpMessage() {
    log.info """
    nf-helminth-amplicon  ${workflow.manifest.version}

    Usage:
      nextflow run . --reference ref.tsv --samplesheet samples.csv --runs runs.tsv
      nextflow run . --reference ref.tsv --samplesheet samples.csv \\
                     --input_dir /data/run01 --platform illumina

    Required:
      --reference        Reference panel TSV (seq_id, locus, type, taxon,
                         parent_id, default_conc, conc_units, corr_factor,
                         sequence)
      --samplesheet      Samplesheet CSV (run_id, platform, sample_id,
                         index_i7, index_i5, sample_type [, file_prefix])

    One of:
      --runs             Run manifest TSV (run_id, platform, path, ...)
      --input_dir        Single-run shortcut; requires --platform

    Options:
      --assay            Assay definition JSON
                         (default: assets/assay/ampseq_helminth.json)
      --index_map        Index primer map (default: assets/index_map.tsv)
      --barcode_repeats  Override the assay's tandem index count (1 or 2)
      --tag_max_edits    Override the assay's tag matching tolerance
      --expected_tag     Override the assay's control tag sequence
      --outdir           Results directory (default: ./results)
    """.stripIndent()
}

workflow {

    if (params.help) {
        helpMessage()
        exit 0
    }

    // Fail on missing inputs here rather than inside a process, so the error
    // arrives immediately and names the parameter.
    def problems = []
    if (!params.reference)   problems << "--reference is required"
    if (!params.samplesheet) problems << "--samplesheet is required"
    if (!params.runs && !params.input_dir) {
        problems << "one of --runs or --input_dir is required"
    }
    if (params.input_dir && !params.platform) {
        problems << "--input_dir requires --platform (illumina or nanopore)"
    }
    if (params.platform && !(params.platform in ['illumina', 'nanopore'])) {
        problems << "--platform must be illumina or nanopore"
    }
    if (problems) {
        log.error "Cannot start:\n  " + problems.join("\n  ")
        exit 1
    }

    ch_reference   = Channel.fromPath(params.reference,   checkIfExists: true)
    ch_samplesheet = Channel.fromPath(params.samplesheet, checkIfExists: true)
    ch_index_map   = Channel.fromPath(params.index_map,   checkIfExists: true)
    ch_assay       = Channel.fromPath(params.assay,       checkIfExists: true)
    ch_runs        = params.runs
        ? Channel.fromPath(params.runs, checkIfExists: true)
        : Channel.fromPath("${projectDir}/assets/NO_FILE", checkIfExists: true)

    ch_primers = params.primers
        ? Channel.fromPath(params.primers, checkIfExists: true)
        : Channel.fromPath("${projectDir}/assets/NO_FILE", checkIfExists: true)
    ch_mocks = params.mocks
        ? Channel.fromPath(params.mocks, checkIfExists: true)
        : Channel.fromPath("${projectDir}/assets/NO_FILE", checkIfExists: true)
    ch_no_file = Channel.fromPath("${projectDir}/assets/NO_FILE", checkIfExists: true)

    PREFLIGHT(ch_reference, ch_samplesheet, ch_runs, ch_index_map, ch_assay,
              ch_primers, ch_mocks)

    // The validation log is worth putting in front of the user rather than
    // leaving in the results directory to be found later.
    PREFLIGHT.out.log
        .map { it.text }
        .view { "\n${it}" }

    // Nothing downstream runs until pre-flight has passed. Everything below
    // depends on its outputs, so Nextflow enforces that ordering without a
    // gate - a failed pre-flight exits non-zero and the DAG stops there.

    // -- what to run on --------------------------------------------------
    // Pre-flight has already resolved every sample to its files and checked
    // they exist, so the samplesheet is read once, here, from its output
    // rather than re-derived. A path that pre-flight blessed is the only
    // path the rest of the pipeline sees.
    ch_resolved = PREFLIGHT.out.file_resolution
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def meta = [
                id        : row.sample_id,
                run_id    : row.run_id,
                platform  : row.platform,
                well      : row.well,
                type      : row.sample_type,
            ]
            [ meta, row.files ? row.files.tokenize(';').collect { file(it) } : [] ]
        }

    ch_illumina = ch_resolved.filter { meta, files -> meta.platform == 'illumina' }
    ch_nanopore_runs = Channel
        .fromPath(params.runs ?: "${projectDir}/assets/NO_FILE")
        .splitCsv(header: true, sep: '\t')
        .filter { row -> row.platform == 'nanopore' }
        .map { row ->
            [ [ run_id: row.run_id, platform: 'nanopore',
                barcode_repeats: row.barcode_repeats ?: null ],
              file(row.path, checkIfExists: true) ]
        }

    // -- the two read paths ----------------------------------------------
    ILLUMINA_PATH(ch_illumina, ch_primers, ch_assay)
    NANOPORE_PATH(ch_nanopore_runs, ch_samplesheet, ch_assay)

    // -- they converge on one assignment ---------------------------------
    ch_seqs = ILLUMINA_PATH.out.seqs.mix(NANOPORE_PATH.out.seqs)

    ASSIGN(
        ch_seqs,
        PREFLIGHT.out.insert_fasta,
        ch_reference.first(),
        ch_primers.first(),
        ch_assay.first(),
    )

    // -- per run ----------------------------------------------------------
    // Synthesis is a run-level step: detection limits come from the run's
    // own blanks and the plate sensors need the whole plate, so the
    // per-sample outputs are collected before it runs.
    ch_run_id = ch_resolved.map { meta, files -> meta.run_id }.first()

    SYNTHESISE(
        ch_run_id,
        ASSIGN.out.counts.collect(),
        ASSIGN.out.summary.collect(),
        ILLUMINA_PATH.out.summaries.collect().ifEmpty([]),
        NANOPORE_PATH.out.summary.ifEmpty(file("${projectDir}/assets/NO_FILE")),
        ch_samplesheet.first(),
        ch_reference.first(),
        PREFLIGHT.out.dir,
    )

    // Amplification efficiency, when the run carried mock communities.
    if (params.mocks) {
        CALIBRATE_MOCKS(
            ch_run_id,
            ch_mocks.first(),
            ASSIGN.out.counts.collect(),
            ch_reference.first(),
        )
    }

    REPORT(
        ch_run_id,
        SYNTHESISE.out.calls
            .mix(SYNTHESISE.out.calibration,
                 SYNTHESISE.out.fates,
                 SYNTHESISE.out.qc)
            .collect(),
        params.comparison
            ? Channel.fromPath(params.comparison, checkIfExists: true)
            : ch_no_file.first(),
    )
}

workflow.onComplete {
    log.info(
        workflow.success
            ? "\nComplete. Report and tables in ${params.outdir}/"
            : "\nFailed. If pre-flight was the failing step, no reads were processed."
    )
}
