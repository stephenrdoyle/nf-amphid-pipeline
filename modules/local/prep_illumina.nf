process PREP_ILLUMINA {
    tag "${meta.id}"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' ?
        'docker://python:3.11-slim' : 'python:3.11-slim' }"

    input:
    tuple val(meta), path(reads)
    path primers
    path assay

    output:
    tuple val(meta), path("prep/${meta.id}.merged.fasta.gz"), emit: merged
    tuple val(meta), path("prep/${meta.id}.pairs.tsv.gz")   , emit: pairs
    tuple val(meta), path("prep/${meta.id}.noprimer.fasta.gz"), emit: noprimer
    path "prep/${meta.id}.prep.json"                        , emit: summary
    path "versions.yml"                                     , emit: versions

    script:
    // Mates arrive as a two-element list in samplesheet order. Handing them
    // over positionally rather than by name keeps this working whatever the
    // sequencing centre called the files.
    def (r1, r2) = [reads[0], reads[1]]
    """
    prep_illumina.py \\
        --r1 ${r1} \\
        --r2 ${r2} \\
        --sample-id ${meta.id} \\
        --primers ${primers} \\
        --assay ${assay} \\
        --min-length ${params.min_length} \\
        --min-overlap ${params.min_overlap} \\
        --max-ee ${params.max_ee} \\
        --outdir prep

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p prep
    echo | gzip > prep/${meta.id}.merged.fasta.gz
    echo | gzip > prep/${meta.id}.pairs.tsv.gz
    echo | gzip > prep/${meta.id}.noprimer.fasta.gz
    echo '{"sample_id":"${meta.id}"}' > prep/${meta.id}.prep.json
    touch versions.yml
    """
}
