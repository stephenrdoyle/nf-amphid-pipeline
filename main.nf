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
      nextflow run . --reference ref.tsv --samplesheet samples.tsv --runs runs.tsv
      nextflow run . --reference ref.tsv --samplesheet samples.tsv \\
                     --input_dir /data/run01 --platform illumina

    Required:
      --reference        Reference panel TSV (seq_id, locus, type, taxon,
                         parent_id, default_conc, conc_units, corr_factor,
                         sequence)
      --samplesheet      Samplesheet CSV (run_id, platform, sample_id,
                         index_i7, index_i5, sample_type [, file_prefix])
      --primers          PCR1 primers TSV (name, direction, sequence)

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
    // Not optional, whatever the placeholder machinery suggests. The
    // primer-trimmed panel is what ASSIGN competes reads against, and
    // pre-flight only builds it from primers - without them ASSIGN never
    // runs, synthesis and the report are skipped with it, and Nextflow
    // prints Complete having analysed nothing.
    if (!params.primers)     problems << "--primers is required (the primer-trimmed panel is built from it)"
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

    // .first() on every one of these, and it is load-bearing, not tidiness.
    //
    // Channel.fromPath on a single file gives a QUEUE channel holding one
    // item. Nextflow pairs a process's inputs positionally and stops at the
    // shortest, so a queue of 8 samples paired with a queue of 1 reference
    // yields ONE task, not 8 - and Nextflow says nothing, because that is
    // the documented behaviour. Seven wells simply never appear, and
    // synthesis reports confidently on the one that did.
    //
    // .first() converts each to a VALUE channel, which is inexhaustible:
    // it is re-read for every sample and imposes no length on the pairing.
    // Doing it here rather than at each call site means a new consumer
    // added later cannot reintroduce the bug.
    ch_reference   = Channel.fromPath(params.reference,   checkIfExists: true).first()
    ch_samplesheet = Channel.fromPath(params.samplesheet, checkIfExists: true).first()
    ch_index_map   = Channel.fromPath(params.index_map,   checkIfExists: true).first()
    ch_assay       = Channel.fromPath(params.assay,       checkIfExists: true).first()
    // Each optional input gets its OWN placeholder. Nextflow stages every
    // input of a process into one directory, so two absent inputs sharing
    // the name NO_FILE collide - "multiple input files for each of the
    // following file names" - and the run dies before the first tool
    // starts. Distinct names cost nothing and make the collision
    // impossible rather than unlikely.
    ch_runs        = ( params.runs
        ? Channel.fromPath(params.runs, checkIfExists: true)
        : Channel.fromPath("${projectDir}/assets/NO_RUNS", checkIfExists: true) ).first()

    ch_primers     = Channel.fromPath(params.primers,     checkIfExists: true).first()
    // Subject metadata is optional and kept OUT of the samplesheet: a
    // samplesheet row is a library and is frozen at prep, a metadata row
    // is a subject and keeps changing, and one subject is several
    // libraries. Pre-flight validates the join rather than assuming it.
    ch_metadata = ( params.metadata
        ? Channel.fromPath(params.metadata, checkIfExists: true)
        : Channel.fromPath("${projectDir}/assets/NO_METADATA", checkIfExists: true) ).first()

    ch_mocks = ( params.mocks
        ? Channel.fromPath(params.mocks, checkIfExists: true)
        : Channel.fromPath("${projectDir}/assets/NO_MOCKS", checkIfExists: true) ).first()
    ch_no_file = Channel.fromPath("${projectDir}/assets/NO_COMPARISON", checkIfExists: true).first()

    PREFLIGHT(ch_reference, ch_samplesheet, ch_runs, ch_index_map, ch_assay,
              ch_primers, ch_mocks, ch_metadata)

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
    // The nanopore runs come from PRE-FLIGHT's resolved manifest, not
    // from the raw runs.tsv.
    //
    // This is the rule stated ten lines above - a path pre-flight
    // blessed is the only path the rest of the pipeline sees - and the
    // nanopore branch used to break it by re-reading the manifest
    // itself. Doing so meant resolving relative paths a second time,
    // in Groovy, and that second implementation was wrong twice: first
    // resolving against the launch directory, then with an isAbsolute()
    // guard that could never fire because Nextflow's file() had already
    // made the path absolute. Both shipped, because Groovy cannot be
    // executed in the environment this is developed in.
    //
    // So the resolution is gone from here rather than fixed again. It
    // happens once, in preflight.py, in Python that has tests.
    ch_nanopore_runs = PREFLIGHT.out.runs_resolved
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
    // A run in which no sample reaches assignment must fail, not
    // succeed quietly. Both branches are conditional - Illumina on the
    // platform column of file_resolution.tsv, nanopore on the platform
    // column of the run manifest - so if those two disagree, or a stub
    // emits rows for neither, every process downstream is simply
    // skipped and Nextflow prints "Complete". That happened: a
    // -stub-run reported success having executed PREFLIGHT and nothing
    // else, and it read as a passing test of the whole DAG.
    ch_seqs = ILLUMINA_PATH.out.seqs
        .mix(NANOPORE_PATH.out.seqs)
        .ifEmpty { error(
            "no sample reached assignment.\n" +
            "  Both platform branches are empty, so nothing after " +
            "pre-flight would run.\n" +
            "  The Illumina branch takes rows of " +
            "preflight/file_resolution.tsv whose platform is " +
            "'illumina';\n" +
            "  the nanopore branch takes rows of the run manifest whose " +
            "platform is 'nanopore'.\n" +
            "  Check that the samplesheet and the run manifest agree on " +
            "the platform column.") }

    ASSIGN(
        ch_seqs,
        // No .first() here, and it is checked rather than assumed.
        // PREFLIGHT is invoked with nothing but value channels, and a
        // process fed only value channels emits value channels - so
        // the panel is already inexhaustible and .first() on it draws
        // a "the operator first is useless" warning. That holds only
        // while every PREFLIGHT input stays a value channel: make one
        // of them a queue and this silently caps ASSIGN at one sample
        // again. tests/test_nextflow.py::test_fanout_is_not_capped
        // works that rule out for itself and fails if it stops being
        // true.
        PREFLIGHT.out.insert_fasta,
        ch_reference,
        ch_primers,
        ch_assay,
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
        NANOPORE_PATH.out.summary.ifEmpty(file("${projectDir}/assets/NO_DEMUX")),
        NANOPORE_PATH.out.split.ifEmpty(file("${projectDir}/assets/NO_SPLIT")),
        ch_samplesheet,
        ch_reference,
        PREFLIGHT.out.dir,
    )

    // Amplification efficiency, when the run carried mock communities.
    if (params.mocks) {
        CALIBRATE_MOCKS(
            ch_run_id,
            ch_mocks,
            ASSIGN.out.counts.collect(),
            ch_reference,
            PREFLIGHT.out.dir,
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
            ? Channel.fromPath(params.comparison, checkIfExists: true).first()
            : ch_no_file,
    )
}

workflow.onComplete {
    log.info(
        workflow.success
            ? "\nComplete. Report and tables in ${params.outdir}/"
            : "\nFailed. If pre-flight was the failing step, no reads were processed."
    )
}
