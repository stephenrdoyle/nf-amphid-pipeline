include { DEMUX_ONT } from '../../modules/local/demux_ont'

workflow NANOPORE_PATH {
    take:
    runs         // [ meta, reads ] - one entry per run, not per sample
    samplesheet
    assay

    main:
    // Demultiplexing is per run: one pooled FASTQ (or MinKNOW directory)
    // in, one file per sample out. Unlike Illumina there is no per-sample
    // parallelism to be had here - the barcodes have to be read before
    // anyone knows which sample a read belongs to.
    DEMUX_ONT(runs, samplesheet, assay)

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
    versions = DEMUX_ONT.out.versions
}
