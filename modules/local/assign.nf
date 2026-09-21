process ASSIGN {
    tag "${meta.id}"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' ?
        'docker://python:3.11-slim' : 'python:3.11-slim' }"

    input:
    tuple val(meta), path(seqs)
    path panel
    path reference
    path primers
    path assay

    output:
    path "assign/${meta.id}.counts.tsv" , emit: counts
    path "assign/${meta.id}.taxa.tsv"   , emit: taxa
    path "assign/${meta.id}.assign.json", emit: summary
    path "versions.yml"                 , emit: versions

    script:
    // Both platforms converge here, but they hand over different things:
    // Illumina a dereplicated FASTA plus a pair table, nanopore a FASTQ
    // whose reads still carry their locus primers. The input shape decides,
    // and the primers are only passed on the path that needs them.
    def input_arg
    if (meta.platform == 'nanopore') {
        input_arg = "--reads ${seqs} --primers ${primers}"
    } else {
        def merged = seqs instanceof List ? seqs.find { it.name.endsWith('.merged.fasta.gz') } : seqs
        def pairs  = seqs instanceof List ? seqs.find { it.name.endsWith('.pairs.tsv.gz') } : null
        input_arg = "--merged ${merged}" + (pairs ? " --pairs ${pairs}" : '')
    }
    """
    assign.py \\
        --panel ${panel} \\
        --reference ${reference} \\
        --sample-id ${meta.id} \\
        ${input_arg} \\
        --assay ${assay} \\
        --max-dist-frac ${params.max_dist_frac} \\
        --min-margin ${params.min_margin} \\
        --outdir assign

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p assign
    touch assign/${meta.id}.counts.tsv assign/${meta.id}.taxa.tsv
    echo '{"sample_id":"${meta.id}"}' > assign/${meta.id}.assign.json
    touch versions.yml
    """
}
