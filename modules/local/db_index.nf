process DB_INDEX {
    tag "${database.name}"
    label 'process_single'
    // Not a result: ~120 MB for a full release, rebuilt in seconds, and
    // useless without the release beside it. Kept out of --outdir.
    publishDir params.outdir, enabled: false

    // No container is declared for the database-mode processes. They need
    // python3 AND minimap2 in one image, and the pipeline does not yet
    // name one; the field laptop runs them natively (docs/USAGE.md).

    input:
    path database

    output:
    // Built once per run and shared by every sample's assignment. The
    // parameters are stored in the index, so assign_db.py cannot map
    // against it with different ones by accident. k=19, w=10: the same
    // accuracy as the map-ont default k=15 on the leave-out benchmark,
    // 1.6x faster; ~1-2 s and ~120 MB for a full release.
    path "refdb.k19w10.mmi", emit: index
    path "versions.yml"    , emit: versions

    script:
    """
    minimap2 -x map-ont -k 19 -w 10 -t ${task.cpus} \\
        -d refdb.k19w10.mmi ${database}/refdb.fasta 2> minimap2_index.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version)
    END_VERSIONS
    """

    stub:
    """
    touch refdb.k19w10.mmi versions.yml
    """
}
