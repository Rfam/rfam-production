nextflow.enable.dsl=2

process tree_files {
    memory '10GB'
    
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
    mysql `python $params.rfamprod/scripts/view/mysql_options.py $params.db` < sql/ftp_rfam_full_reqion.sql > Rfam.full_region
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
    mysql `python $params.rfamprod/scripts/view/mysql_options.py $params.db` < sql/ftp_rfam_pdb.sql > Rfam.pdb
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

process generate_fasta_files {
    memory '10GB'
    publishDir "$params.release_ftp", mode: "copy"
   
    output:
    val('done')

    """
    mkdir $params.release_ftp/fasta_files
    python $params.rfamprod/scripts/export/fasta_file_generator.py --seq-db $params.rfamseqfa --rfam-seed $params.release_ftp/seed/Rfam.seed --all --outdir $params.release_ftp/fasta_files
    """
}


workflow generate_ftp_files {
    tree_files()
    generate_full_region_file \
    | generate_pdb_file \
    | generate_clanin_file
    generate_fasta_files()
}

workflow {
    generate_ftp_files()
}
