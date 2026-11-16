process LIMA_PACBIO {
    tag "${meta.run_id}"
    label 'process_medium'

    // lima (PacBio's demultiplexer) and python3 together; no image names
    // both yet, so as for the database-mode processes none is declared.
    // lima has no macOS build: main.nf refuses --pacbio_demux lima on a
    // Mac, and this refuses again where the task actually runs.
    conda "bioconda::lima=2.9.0 conda-forge::python=3.11"

    input:
    tuple val(meta), path(reads)
    path samplesheet
    path primers
    path assay

    output:
    tuple val(meta), path("demux/reads/*.fastq.gz"), emit: reads
    path "demux/demux_summary.json"                , emit: summary
    path "demux/index_matrix.tsv"                  , emit: matrix
    path "versions.yml"                            , emit: versions

    script:
    // lima instead of the built-in demultiplexer: the run's tags written
    // as lima barcodes with the preset that fits them (demux_pacbio.py
    // lima_design: one barcode per sample and SYMMETRIC where each
    // sample's two tags are one or differ only at the primer-adjacent
    // base; one per tag and ASYMMETRIC where none are; a mix refused),
    // one lima call per input file, and its per-pair outputs
    // read back as samples by demux_pacbio.py --from-lima, which writes
    // the same summary and matrix as the built-in route.
    """
    if [ "\$(uname -s)" = "Darwin" ]; then
        echo "lima has no macOS build: run on Linux (the farm), or use --pacbio_demux builtin" >&2
        exit 1
    fi
    demux_pacbio.py \\
        --samplesheet ${samplesheet} \\
        --run-id ${meta.run_id} \\
        --primers ${primers} \\
        --write-lima-barcodes barcodes.fasta

    if [ -d ${reads} ]; then
        inputs=\$(find -L ${reads} -maxdepth 1 -type f \\( -name '*.bam' -o -name '*.fastq*' -o -name '*.fq*' \\) | sort)
    else
        inputs=${reads}
    fi
    mkdir -p lima
    k=0
    for f in \$inputs; do
        k=\$((k + 1))
        case "\$f" in *.bam) ext=bam ;; *) ext=fastq ;; esac
        lima "\$f" barcodes.fasta lima/part\$k.\$ext \\
            --hifi-preset \$(cat barcodes.fasta.preset) \\
            --split-named \\
            --num-threads ${task.cpus}
    done

    demux_pacbio.py \\
        --samplesheet ${samplesheet} \\
        --run-id ${meta.run_id} \\
        --primers ${primers} \\
        --from-lima lima \\
        --outdir demux

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lima: \$(lima --version | head -1 | sed 's/lima //')
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
