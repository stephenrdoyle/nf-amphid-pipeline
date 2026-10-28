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
include { DB_INDEX       } from './modules/local/db_index'
include { ASSIGN_DB      } from './modules/local/assign_db'
include { SYNTHESISE     } from './modules/local/synthesise'
include { CALIBRATE_MOCKS} from './modules/local/calibrate_mocks'
include { COMPARE_RUNS   } from './modules/local/compare_runs'
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
        or
      --reference_database
                         A vespa-18s-refdb release directory: classify
                         against a taxonomic database instead of a panel
                         (nanopore only so far; needs minimap2)
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
      --spiked           The quantification standards were added to every
                         sample: calibrate, give copy numbers, and fail
                         wells missing them (default: unspiked)
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
    // A panel or a database release: one of them, never both. The
    // database route classifies against a release rather than competing
    // reads across a panel (docs/REFERENCE_DATABASE.md).
    if (!params.reference && !params.reference_database) {
        problems << "one of --reference (a panel) or --reference_database (a release) is required"
    }
    if (params.reference && params.reference_database) {
        problems << "--reference and --reference_database are alternatives: give one"
    }
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
    // A flag. Compared as text, never by truthiness: from Nextflow 26.04
    // every command-line parameter arrives as a string, so a bare
    // --spiked may be 'true' rather than true, and --spiked false is the
    // string 'false' - truthy in Groovy, so it would quietly mean spiked.
    // Anything but true or false (e.g. --spiked yes) is refused.
    if (!(params.spiked.toString().toLowerCase() in ['true', 'false'])) {
        problems << "--spiked is a flag: give it with no value (got '${params.spiked}')"
    }
    // The test profiles run on simulated reads, which are generated rather
    // than committed - so a fresh clone has none. Without this check the
    // first run on a new machine dies inside pre-flight on "path does not
    // exist", and then the empty-branch guard adds a second, misleading
    // error about platforms. Each entry is the reads directory a profile's
    // runs.tsv points at, and the command that creates it.
    def fixtures = [
        'test'    : [ "${projectDir}/tests/data/sim/RUN01",
                      "python3 bin/simulate.py --reference tests/data/reference.tsv " +
                      "--samplesheet tests/data/samplesheet.tsv --runs tests/data/runs.tsv " +
                      "--assay assets/assay/ampseq_helminth.json --outdir tests/data/sim " +
                      "--reads-per-sample 300" ],
        'test_ont': [ "${projectDir}/tests/data/ont/run",
                      "python3 tests/make_ont_fixture.py" ],
        'test_db' : [ "${projectDir}/tests/data/db/run",
                      "python3 tests/make_db_fixture.py" ],
    ]
    (workflow.profile ?: '').tokenize(',').each { prof ->
        def fx = fixtures[prof.trim()]
        if (fx) {
            def reads = file(fx[0])
            if (!reads.isDirectory() || !reads.list()) {
                problems << "-profile ${prof.trim()}: the fixture reads have not been " +
                            "generated (${fx[0]} is missing or empty).\n" +
                            "  They are simulated, not committed. From ${projectDir}, run:\n" +
                            "    ${fx[1]}"
            }
        }
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
    // In database mode there is no panel: pre-flight gets the
    // NO_REFERENCE placeholder in its place, and the release instead of
    // NO_DATABASE. Each has its own placeholder name (see below).
    ch_reference   = ( params.reference
        ? Channel.fromPath(params.reference, checkIfExists: true)
        : Channel.fromPath("${projectDir}/assets/NO_REFERENCE", checkIfExists: true) ).first()
    ch_database    = ( params.reference_database
        ? Channel.fromPath(params.reference_database, type: 'dir', checkIfExists: true)
        : Channel.fromPath("${projectDir}/assets/NO_DATABASE", checkIfExists: true) ).first()
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
              ch_primers, ch_mocks, ch_metadata, ch_database)

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
            // id is the SANITISED form, not the declared one. It is
            // interpolated unquoted into shell commands (--sample-id) and
            // names every per-sample file, so a declared 'S 1' split the
            // argument and 'A/B' made a directory. Pre-flight has already
            // refused ids whose sanitised forms collide, and synthesis
            // restores the declared id through file_resolution.tsv - the
            // same route the nanopore branch's ids have always taken.
            def meta = [
                id        : row.sanitised_id,
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
                barcode_repeats: row.barcode_repeats ?: null,
                index_length: row.index_length ?: null ],
              file(row.path, checkIfExists: true) ]
        }

    // -- the two read paths ----------------------------------------------
    ILLUMINA_PATH(ch_illumina, ch_primers, ch_assay)
    NANOPORE_PATH(ch_nanopore_runs, ch_samplesheet, ch_assay,
                  PREFLIGHT.out.lookalikes)

    // -- they converge on one assignment ---------------------------------
    // A run in which no sample reaches assignment must fail, not
    // succeed quietly. Both branches are conditional - Illumina on the
    // platform column of file_resolution.tsv, nanopore on the platform
    // column of the run manifest - so if those two disagree, or a stub
    // emits rows for neither, every process downstream is simply
    // skipped and Nextflow prints "Complete". That happened: a
    // -stub-run reported success having executed PREFLIGHT and nothing
    // else, and it read as a passing test of the whole DAG.
    //
    // The guard fires only once pre-flight has PASSED, and not on a run
    // that is being stopped. As a bare ifEmpty it also fired whenever
    // pre-flight failed or the run was interrupted - every branch is
    // empty then too - and added a second, wrong explanation about
    // platform columns under the real one.
    ch_seqs = ILLUMINA_PATH.out.seqs
        .mix(NANOPORE_PATH.out.seqs)
    PREFLIGHT.out.summary
        .combine(ch_seqs.count())
        .subscribe { _summary, n ->
            if (n == 0 && !workflow.session.isAborted()) {
                error(
                    "no sample reached assignment.\n" +
                    "  Both platform branches are empty, so nothing after " +
                    "pre-flight would run.\n" +
                    "  The Illumina branch takes rows of " +
                    "preflight/file_resolution.tsv whose platform is " +
                    "'illumina';\n" +
                    "  the nanopore branch takes rows of the run manifest " +
                    "whose platform is 'nanopore'.\n" +
                    "  Check that the samplesheet and the run manifest agree " +
                    "on the platform column.")
            }
        }

    if (params.reference_database) {
        // -- database mode: classify consensus sequences ------------------
        // One minimap2 index for the run, then one assignment per sample.
        // The index is built from a value channel only, so its output is
        // a value channel and pairs with every sample (see the .first()
        // note above).
        DB_INDEX(ch_database)
        ASSIGN_DB(
            ch_seqs,
            ch_database,
            DB_INDEX.out.index,
            ch_primers,
        )
        // Synthesis and the report for database mode are not built yet:
        // the per-sample cluster and taxon tables in assign_db/ are the
        // result. Said at the end of the run, so nobody waits for a
        // report that will not come.
        ASSIGN_DB.out.summary.collect().subscribe { files ->
            log.info "\nDatabase mode: ${files.size()} sample(s) assigned. " +
                     "Per-sample clusters and taxa are in ${params.outdir}/assign_db/; " +
                     "synthesis and the report are not built for database mode yet."
        }
    } else {
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

        // Cross-run comparison. --comparison is a previous run's calls.tsv,
        // and report.py reads the comparison.json that compare_runs.py makes
        // from it - so it is compared here, not handed to the report raw.
        // Passing the TSV straight through, as this once did, gave the report
        // a file it could not parse. The previous run's label is left empty
        // so compare_runs.py takes it from that file's own run_id column.
        ch_comparison = ch_no_file
        if (params.comparison) {
            COMPARE_RUNS(
                SYNTHESISE.out.calls
                    .combine(ch_run_id)
                    .combine(Channel.fromPath(params.comparison, checkIfExists: true))
                    .map { calls, rid, previous -> [ rid, calls, '', previous ] }
            )
            ch_comparison = COMPARE_RUNS.out.json.first()
        }

        // Every synthesis output report.py reads must be staged here - the
        // report renders a section only when its file exists, so a missing
        // one drops the section without a word. group_summary.tsv was missed
        // that way. It is optional, and mix() of an empty channel adds
        // nothing, so a run with no grouping still reports.
        // tests/test_nextflow.py::test_report_gets_what_it_reads holds this.
        REPORT(
            ch_run_id,
            SYNTHESISE.out.calls
                .mix(SYNTHESISE.out.calibration,
                     SYNTHESISE.out.fates,
                     SYNTHESISE.out.qc,
                     SYNTHESISE.out.groups)
                .collect(),
            ch_comparison,
        )
    }
}

workflow.onComplete {
    log.info(
        workflow.success
            ? "\nComplete. Report and tables in ${params.outdir}/"
            : "\nFailed. If pre-flight was the failing step, no reads were processed."
    )
}
