include { PREP_ILLUMINA } from '../../modules/local/prep_illumina'

workflow ILLUMINA_PATH {
    take:
    samples      // [ meta, [r1, r2] ]
    primers
    assay

    main:
    // Per sample: the sequencer already demultiplexed, so there is nothing
    // run-wide to do and every well is independent.
    PREP_ILLUMINA(samples, primers, assay)

    // Assignment wants the merged reads and the unmerged pairs together -
    // a pair that could not be merged is still evidence, assigned with
    // both mates required to agree.
    seqs = PREP_ILLUMINA.out.merged
        .join(PREP_ILLUMINA.out.pairs)
        .map { meta, merged, pairs -> [ meta, [ merged, pairs ] ] }

    emit:
    seqs     = seqs
    summaries = PREP_ILLUMINA.out.summary
    versions = PREP_ILLUMINA.out.versions
}
