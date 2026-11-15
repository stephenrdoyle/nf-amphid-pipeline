include { DEMUX_PACBIO } from '../../modules/local/demux_pacbio'
include { LIMA_PACBIO  } from '../../modules/local/lima_pacbio'

workflow PACBIO_PATH {
    take:
    runs         // [ meta, reads ] - one entry per run (one library)
    samplesheet
    primers
    assay

    main:
    // Demultiplexed once per run, by the run's own tags: a SMRT cell
    // shared between projects carries other libraries' tags, and those
    // reads must not reach this run's samples. The built-in
    // demultiplexer by default; lima (Linux only) when asked.
    if (params.pacbio_demux == 'lima') {
        LIMA_PACBIO(runs, samplesheet, primers, assay)
        ch_reads    = LIMA_PACBIO.out.reads
        ch_summary  = LIMA_PACBIO.out.summary
        ch_matrix   = LIMA_PACBIO.out.matrix
        ch_versions = LIMA_PACBIO.out.versions
    } else {
        DEMUX_PACBIO(runs, samplesheet, primers, assay)
        ch_reads    = DEMUX_PACBIO.out.reads
        ch_summary  = DEMUX_PACBIO.out.summary
        ch_matrix   = DEMUX_PACBIO.out.matrix
        ch_versions = DEMUX_PACBIO.out.versions
    }

    // One channel item per sample, the id the whole filename less its
    // suffix (not simpleName, which stops at the first dot), as for
    // nanopore.
    seqs = ch_reads
        .transpose()
        .map { meta, fastq ->
            [ meta + [ id: fastq.name.replaceFirst(/\.fastq\.gz$/, ''),
                       platform: 'pacbio' ], fastq ]
        }

    emit:
    seqs     = seqs
    summary  = ch_summary
    matrix   = ch_matrix
    versions = ch_versions
}
