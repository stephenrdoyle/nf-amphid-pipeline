process SPLIT_CONCATEMERS {
    tag "${meta.run_id}:${index}"
    label 'process_shard'

    container "${ workflow.containerEngine == 'singularity' ?
        'docker://python:3.11-slim' : 'python:3.11-slim' }"

    input:
    tuple val(meta), val(index), path(reads)
    path assay
    path lookalikes

    output:
    // One output FASTQ per input FASTQ, which is why there is no merge
    // afterwards: demultiplexing walks a directory, so the shards ARE
    // the output. Collecting them is staging, not a rewrite of tens of
    // gigabytes.
    tuple val(meta), path("split/reads/*.fastq.gz"), emit: reads
    path "split/parts/*.json"                      , emit: part
    path "versions.yml"                            , emit: versions

    script:
    def repeats = meta.barcode_repeats ? "--barcode-repeats ${meta.barcode_repeats}" : ''
    // The index span sets where each cut falls. Pre-flight reads it off
    // the run's declared indices; left to the tool's default of 8, a 10bp
    // index set put every cut 4bp off.
    def idx_len = meta.index_length ? "--index-length ${meta.index_length}" : ''
    // One input file per task, so no --threads: the tool divides work
    // BY INPUT FILE, and a shard holding one file has nothing for a
    // second worker to do. Asking for more cores would reserve them
    // and leave them idle - on a 74 file run, over a thousand cores
    // held for nothing, in slots LSF finds far harder to place than
    // the single cores this actually needs.
    """
    split_concatemers.py \\
        --reads ${reads} \\
        --assay ${assay} \\
        --part-index ${index} \\
        --lookalikes ${lookalikes} \\
        ${repeats} ${idx_len} \\
        --anchor-max-edits ${params.anchor_max_edits} \\
        --min-length ${params.min_unit_length} \\
        --outdir split

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p split/reads split/parts
    echo | gzip > split/reads/${index}_stub.fastq.gz
    echo '{"source":"stub_${index}","reads_in":0,"pieces_out":0,"reads_split":0,"pieces_dropped_below_min_length":0,"pieces_per_read_distribution":{},"bases_in":0,"bases_out":0}' > split/parts/${index}_stub.json
    touch versions.yml
    """
}
