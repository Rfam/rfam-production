nextflow.enable.dsl=2

process tree_files {
    memory '10GB'
    input:
    val(_flag)
    
    output:
    val('tree_done')

    """
    rm -rf $params.release_ftp/tree/
    mkdir -p $params.release_ftp/tree/
    python $params.rfamprod/scripts/export/generate_ftp_files.py --acc all --tree --dest-dir $params.release_ftp/tree
    """
}

process generate_full_region_file {
    publishDir "$params.release_ftp", mode: "copy"
       
    output:
    path('Rfam.full_region.gz')

    """
    mysql `python $params.rfamprod/scripts/view/mysql_options.py $params.db` < $params.rfamprod/sql/ftp_rfam_full_reqion.sql > Rfam.full_region
    gzip Rfam.full_region
    """

}

process generate_pdb_file {
    publishDir "$params.release_ftp", mode: "copy"

    input:
    path(query)

    output:
    path('Rfam.pdb.gz')

    """
    mysql `python $params.rfamprod/scripts/view/mysql_options.py $params.db` < $params.rfamprod/sql/ftp_rfam_pdb.sql > Rfam.pdb
    gzip Rfam.pdb
    """

}

process generate_clanin_file {
    publishDir "$params.release_ftp", mode: "copy"

    input:
    path(query)
    
    output:
    val('done')

    """
    python $params.rfamprod/scripts/release/clanin_file_generator.py $params.release_ftp
    """
}

process generate_3d_seed_file {

    input:
    path(query)

    output:
    val('done')

    """
    python $params.rfamprod/scripts/release/generate_3d_seed_file.py -f $params.3d_families --seed $params.release_ftp/seed
    """
}


workflow generate_ftp_files {
    take: start 
    emit: done
    main:
        tree_files(start)
        generate_full_region_file \
        | generate_pdb_file \
        | generate_clanin_file \
        | generate_3d_seed_file
        | set { done }
}

workflow {
    generate_ftp_files(Channel.of('start'))
}
