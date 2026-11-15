process DEMUX_PACBIO {
    tag "${meta.run_id}"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' ?
        'docker://python:3.11-slim' : 'python:3.11-slim' }"

    input:
    tuple val(meta), path(reads)
    path samplesheet
    path primers
    path assay

    output:
    // One FASTQ per sample with reads, tags off and primers on; the
    // protocol's own trim cuts the primers (docs/CONTRACTS.md).
    tuple val(meta), path("demux/reads/*.fastq.gz"), emit: reads
    path "demux/demux_summary.json"                , emit: summary
    path "demux/index_matrix.tsv"                  , emit: matrix
    path "versions.yml"                            , emit: versions

    script:
    // The built-in demultiplexer: the library's FASTQ or BAM (a file, or
    // a directory of a library's SMRT cells), matched to the run's tags
    // from the samplesheet (index_i5 forward, index_i7 reverse). The
    // assay JSON's optional "pacbio" section sets the tolerances.
    def tag_arg = params.pacbio_tag_max_edits != null ? "--tag-max-edits ${params.pacbio_tag_max_edits}" : ''
    """
    demux_pacbio.py \\
        --reads ${reads} \\
        --samplesheet ${samplesheet} \\
        --run-id ${meta.run_id} \\
        --primers ${primers} \\
        --assay ${assay} \\
        ${tag_arg} \\
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
    touch demux/index_matrix.tsv
    echo '{"run_id":"${meta.run_id}"}' > demux/demux_summary.json
    touch versions.yml
    """
}
