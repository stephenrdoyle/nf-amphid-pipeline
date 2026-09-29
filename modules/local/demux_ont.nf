process DEMUX_ONT {
    tag "${meta.run_id}"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' ?
        'docker://python:3.11-slim' : 'python:3.11-slim' }"

    input:
    tuple val(meta), path(reads)
    path samplesheet
    path assay

    output:
    // One FASTQ per sample. The run is demultiplexed once, so this is a
    // single process emitting many files rather than one task per sample.
    tuple val(meta), path("demux/reads/*.fastq.gz"), emit: reads
    path "demux/assignments.tsv"                   , emit: assignments
    path "demux/demux_summary.json"                , emit: summary
    path "demux/index_matrix.tsv"                  , emit: matrix, optional: true
    path "versions.yml"                            , emit: versions

    script:
    // A MinION run is a directory of batched FASTQ, not a file. The tool
    // walks it, uses fastq_pass and reports what it excluded, so whatever
    // was staged - a file, several files, or a run folder - is passed
    // through unchanged.
    def repeats = meta.barcode_repeats ? "--barcode-repeats ${meta.barcode_repeats}" : ''
    def failed  = params.include_failed ? '--include-failed' : ''
    """
    demux_ont.py \\
        --reads ${reads} \\
        --samplesheet ${samplesheet} \\
        --run-id ${meta.run_id} \\
        --assay ${assay} \\
        ${repeats} ${failed} \\
        --max-mismatch ${params.max_mismatch} \\
        --threads ${task.cpus} \\
        --outdir demux

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p demux/reads
    echo | gzip > demux/reads/stub.fastq.gz
    touch demux/assignments.tsv demux/index_matrix.tsv
    echo '{"run_id":"${meta.run_id}"}' > demux/demux_summary.json
    touch versions.yml
    """
}
