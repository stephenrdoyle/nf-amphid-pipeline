include { SPLIT_LIST        } from '../../modules/local/split_list'
include { SPLIT_CONCATEMERS } from '../../modules/local/split_concatemers'
include { MERGE_SPLIT       } from '../../modules/local/merge_split'
include { DEMUX_ONT         } from '../../modules/local/demux_ont'

workflow NANOPORE_PATH {
    take:
    runs         // [ meta, reads ] - one entry per run, not per sample
    samplesheet
    assay

    main:
    // Cut first. Ligation joins amplicons into single reads, and an
    // uncut join presents one molecule's i5 with another's i7 - a pair
    // of real indices belonging to a well that held neither. Once that
    // pair has been looked up there is no evidence left that anything
    // was wrong, so the cut has to happen before demultiplexing, not
    // after.
    //
    // Splitting is the stage that walks every base of every read, and
    // on a full flow cell it is hours of work - it is the step that
    // hit a twelve hour wall clock limit on the first real run. So it
    // is fanned out here rather than threaded inside one task: a
    // MinKNOW run arrives as dozens of batch files, the work divides
    // by file with nothing shared between them, and one single-core
    // task per file lets LSF pack them into whatever is free. The
    // alternative - one sixteen-core task - is bounded by one node and
    // waits for a sixteen-core slot to open.
    SPLIT_LIST(runs, assay)

    // The tool decided which files the run consists of; this only
    // turns its answer into one channel item per file, carrying the
    // index so each shard names its output as the whole-run route
    // would have. Without it every shard numbers its file 0000 and two
    // MinKNOW batches sharing a basename collide on collection.
    ch_shards = SPLIT_LIST.out.list
        .splitCsv(header: false, sep: '\t', elem: 1)
        .map { meta, row -> [ meta, row[0], file(row[1], checkIfExists: true) ] }

    SPLIT_CONCATEMERS(ch_shards, assay)

    // No merge of the reads themselves: the splitter writes one FASTQ
    // per input FASTQ and demultiplexing reads a directory, so
    // gathering the shards is Nextflow staging them side by side. The
    // meta is the run, identical across shards, so groupTuple collects
    // the run's pieces back into one item.
    ch_split_reads = SPLIT_CONCATEMERS.out.reads
        .groupTuple()
        .map { meta, files -> [ meta, files.flatten() ] }

    // The summary IS merged, because every read-funnel percentage
    // downstream is taken against its reads_in.
    MERGE_SPLIT(
        SPLIT_LIST.out.list.map { meta, f -> meta }.first(),
        SPLIT_CONCATEMERS.out.part.collect(),
        SPLIT_LIST.out.report,
        ch_shards.count(),
    )

    // Demultiplexing is per run: the pooled reads in, one file per
    // sample out. Unlike splitting there is no per-file parallelism to
    // be had - a sample's reads are scattered across every batch file,
    // so the gather would be keyed by sample rather than by file.
    DEMUX_ONT(ch_split_reads, samplesheet, assay)

    // Fan out to one channel item per sample, carrying the run's metadata
    // forward so assignment knows it is on the nanopore path.
    seqs = DEMUX_ONT.out.reads
        .transpose()
        .map { meta, fastq ->
            [ meta + [ id: fastq.simpleName, platform: 'nanopore' ], fastq ]
        }

    emit:
    seqs     = seqs
    summary  = DEMUX_ONT.out.summary
    // Carried through to synthesis so the read funnel can start at reads
    // off the sequencer. Without it the funnel starts at amplicon units,
    // and every percentage below is taken against a denominator the
    // split quietly changed.
    split    = MERGE_SPLIT.out.summary
    versions = DEMUX_ONT.out.versions
        .mix(SPLIT_CONCATEMERS.out.versions,
             SPLIT_LIST.out.versions,
             MERGE_SPLIT.out.versions)
}
