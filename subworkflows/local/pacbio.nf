include { DEMUX_PACBIO } from '../../modules/local/demux_pacbio'

workflow PACBIO_PATH {
    take:
    runs         // [ meta, reads ] - one entry per run (one library)
    samplesheet
    primers
    assay

    main:
    // Demultiplexed once per run, by the run's own tags: a SMRT cell
    // shared between projects carries other libraries' tags, and those
    // reads must not reach this run's samples.
    DEMUX_PACBIO(runs, samplesheet, primers, assay)

    // One channel item per sample, the id the whole filename less its
    // suffix (not simpleName, which stops at the first dot), as for
    // nanopore.
    seqs = DEMUX_PACBIO.out.reads
        .transpose()
        .map { meta, fastq ->
            [ meta + [ id: fastq.name.replaceFirst(/\.fastq\.gz$/, ''),
                       platform: 'pacbio' ], fastq ]
        }

    emit:
    seqs     = seqs
    summary  = DEMUX_PACBIO.out.summary
    matrix   = DEMUX_PACBIO.out.matrix
    versions = DEMUX_PACBIO.out.versions
}
