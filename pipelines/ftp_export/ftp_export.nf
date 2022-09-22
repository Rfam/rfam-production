process tree_files {
    memory '10GB'
    input:
    val(_flag)

    output:
    val('tree_done')

    """
    rm -rf $params.ftp_exports/tree/
    mkdir -p $params.ftp_exports/tree/
    python $params.scripts/export/generate_ftp_files.py --acc all --tree --dest-dir $params.ftp_exports/tree
    """
}

process generate_full_region_file {
    publishDir "$params.ftp_exports", mode: "copy"

    output:
    path('Rfam.full_region.gz')

    """
    mysql `python $params.scripts/view/mysql_options.py $params.db` < $params.sql/ftp_rfam_full_reqion.sql > Rfam.full_region
    gzip Rfam.full_region
    """

}

process generate_pdb_file {
    publishDir "$params.ftp_exports", mode: "copy"

    input:
    path(query)

    output:
    path('Rfam.pdb.gz')

    """
    mysql `python $params.scripts/view/mysql_options.py $params.db` < $params.sql/ftp_rfam_pdb.sql > Rfam.pdb
    gzip Rfam.pdb
    """

}

process generate_clanin_file {
    publishDir "$params.ftp_exports", mode: "copy"

    input:
    path(query)

    output:
    val('done')

    """
    python $params.scripts/release/clanin_file_generator.py $params.ftp_exports
    """
}

process generate_fasta_files {
    memory '10GB'

    output:
    val('done')

    """
    rm -rf mkdir $params.ftp_exports/fasta_files
    mkdir $params.ftp_exports/fasta_files
    python $params.scripts/export/fasta_file_generator.py --seq-db $params.rfamseqfa --rfam-seed $params.ftp_exports/seed/Rfam.seed --all --outdir $params.ftp_exports/fasta_files
    """
}


workflow generate_ftp_files {
    take: start
    emit: done
    main:
        tree_files(start)
        generate_full_region_file \
        | generate_pdb_file \
        | generate_clanin_file
        generate_fasta_files() \
        | set { done }
}

workflow {
    generate_ftp_files(Channel.of('start'))
}
