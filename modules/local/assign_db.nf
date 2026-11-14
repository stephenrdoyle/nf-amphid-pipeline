process ASSIGN_DB {
    tag "${meta.id}"
    label 'process_db'

    // See db_index.nf: no container until one image carries python3 and
    // minimap2 together.

    input:
    tuple val(meta), path(reads)
    path database
    path index
    path primers

    output:
    path "assign_db/${meta.id}.clusters.tsv"   , emit: clusters
    path "assign_db/${meta.id}.taxa.tsv"       , emit: taxa
    path "assign_db/${meta.id}.consensus.fasta", emit: consensus
    path "assign_db/${meta.id}.assign_db.json" , emit: summary
    path "versions.yml"                        , emit: versions

    script:
    // Per sample, single-threaded under the field profile like ASSIGN:
    // the samples are the parallelism. minimap2 gets the task's cores
    // where a profile gives it more than one.
    """
    assign_db.py \\
        --reads ${reads} \\
        --database ${database} \\
        --index ${index} \\
        --primers ${primers} \\
        --sample-id ${meta.id} \\
        --threads ${task.cpus} \\
        --min-cluster-reads ${params.db_min_cluster_reads} \\
        --lca-margin ${params.db_lca_margin} \
        --primer-max-edits ${params.db_primer_max_edits} \\
        --floors '${params.db_floors}' \\
        --calling-rule ${params.db_calling_rule} \\
        --gap-floors '${params.db_gap_floors}' \\
        --gap-edits ${params.db_gap_edits} \\
        --min-species-span ${params.db_min_species_span} \\
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
