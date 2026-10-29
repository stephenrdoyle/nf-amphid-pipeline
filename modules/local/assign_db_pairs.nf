process ASSIGN_DB_PAIRS {
    tag "${meta.id}"
    label 'process_low'

    // See db_index.nf: no container until one image carries python3 and
    // minimap2 together.

    input:
    tuple val(meta), path(reads)
    path database
    path index
    path primers
    path assay

    output:
    path "assign_db/${meta.id}.clusters.tsv"   , emit: clusters
    path "assign_db/${meta.id}.taxa.tsv"       , emit: taxa
    path "assign_db/${meta.id}.consensus.fasta", emit: consensus
    path "assign_db/${meta.id}.assign_db.json" , emit: summary
    path "versions.yml"                        , emit: versions

    script:
    // Database mode on Illumina: the pairs are not merged (at 2x150 no V4
    // insert overlaps), but windowed, denoised and classified as pairs.
    // Outputs have assign_db.py's shape, so synthesis reads either
    // platform. Windows default to the release's illumina_windows.
    def (r1, r2) = [reads[0], reads[1]]
    def w1 = params.db_window_r1 != null ? "--window-r1 ${params.db_window_r1}" : ''
    def w2 = params.db_window_r2 != null ? "--window-r2 ${params.db_window_r2}" : ''
    """
    assign_db_pairs.py \\
        --r1 ${r1} \\
        --r2 ${r2} \\
        --database ${database} \\
        --index ${index} \\
        --primers ${primers} \\
        --assay ${assay} \\
        --sample-id ${meta.id} \\
        --threads ${task.cpus} \\
        ${w1} ${w2} \
        --read-length ${params.illumina_read_length} \\
        --max-ee ${params.max_ee} \\
        --lca-margin ${params.db_lca_margin} \
        --primer-max-edits ${params.db_primer_max_edits} \\
        --floors '${params.db_floors_illumina}' \\
        --outdir assign_db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        minimap2: \$(minimap2 --version)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p assign_db
    touch assign_db/${meta.id}.clusters.tsv assign_db/${meta.id}.taxa.tsv
    touch assign_db/${meta.id}.consensus.fasta
    echo '{"sample_id":"${meta.id}"}' > assign_db/${meta.id}.assign_db.json
    touch versions.yml
    """
}
