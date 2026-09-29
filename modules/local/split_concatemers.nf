process SPLIT_CONCATEMERS {
    tag "${meta.run_id}"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' ?
        'docker://python:3.11-slim' : 'python:3.11-slim' }"

    input:
    tuple val(meta), path(reads)
    path assay

    output:
    // One FASTQ of amplicon units, replacing the pooled reads. Runs
    // before demultiplexing rather than after, because the whole point
    // is that the index pair on an uncut read may belong to two
    // different molecules - there is nothing to fix once it has been
    // demultiplexed on that pair.
    tuple val(meta), path("split/reads/*.fastq.gz"), emit: reads
    path "split/split_summary.json"              , emit: summary
    path "versions.yml"                          , emit: versions

    script:
    def repeats = meta.barcode_repeats ? "--barcode-repeats ${meta.barcode_repeats}" : ''
    """
    split_concatemers.py \\
        --reads ${reads} \\
        --assay ${assay} \\
        ${repeats} \\
        --anchor-max-edits ${params.anchor_max_edits} \\
        --min-length ${params.min_unit_length} \\
        --threads ${task.cpus} \\
        --outdir split

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p split/reads
    echo | gzip > split/reads/stub.fastq.gz
    echo '{"reads_in":0,"pieces_out":0}' > split/split_summary.json
    touch versions.yml
    """
}
