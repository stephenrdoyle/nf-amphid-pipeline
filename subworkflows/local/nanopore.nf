include { SPLIT_CONCATEMERS } from '../../modules/local/split_concatemers'
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
    SPLIT_CONCATEMERS(runs, assay)

    // Demultiplexing is per run: one pooled FASTQ (or MinKNOW directory)
    // in, one file per sample out. Unlike Illumina there is no per-sample
    // parallelism to be had here - the barcodes have to be read before
    // anyone knows which sample a read belongs to.
    DEMUX_ONT(SPLIT_CONCATEMERS.out.reads, samplesheet, assay)

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
    split    = SPLIT_CONCATEMERS.out.summary
    versions = DEMUX_ONT.out.versions.mix(SPLIT_CONCATEMERS.out.versions)
}
