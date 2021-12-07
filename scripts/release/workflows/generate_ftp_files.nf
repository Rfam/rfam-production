nextflow.enable.dsl=2

process tree_files {
    memory '10GB'
    
    output:
    val('tree_done')

    """
    rm -rf ${params.release_ftp}/tree/
    mkdir -p ${params.release_ftp}/tree/
    python ${params.rfamprod}/scripts/export/generate_ftp_files.py --acc all --seed --dest-dir ${params.release_ftp}/tree
    """
}

process rfam2go_file {
    input:
    val('tree_done')
    
    output:
    val('rfam2go_export')

    """
    rm -rf ${params.release_ftp}/rfam2go
    mkdir -p ${params.release_ftp}/rfam2go
    perl ${params.perl_path}/export/rfam2go.pl > ${params.release_ftp}/rfam2go
    """
}

process add_header {
    input:
    val('rfam2go_export')
    
    output:
    val('rfam2go_header')

    """
    python ${params.rfamprod}/scripts/release/add_header_rfam2go.py --input $query
    """
}

process checksum {
    publishDir "${params.release_ftp}/rfam2go/", mode: "copy"
    
    input:
    val('rfam2go_header')
    
    output:
    path("md5.txt")

    """
    cd ${params.release_ftp}/rfam2go && md5sum * > md5.txt
    """
}

process validate_go {
    input:
    path(query)
    
    output:
    val('go_validated')

    """
    perl ${params.perl_path}/qc/validate_rfam2go.pl goselect selectgo ${params.release_ftp}/rfam2go
    """

}

process generate_full_region_file {
    publishDir "${params.release_ftp}", mode: "copy"
    input:
    val('go_validated')
    
    output:
    path("Rfam.full_region.gz")


    """
    mysql `python $params.rfamprod/scripts/view/mysql_options.py $params.db` < sql/ftp_rfam_full_reqion.sql > Rfam.full_region
    gzip Rfam.full_region
    """

}

process generate_pdb_file {
    input:
    path("Rfam.full_region.gz")

    
    output:
    path("Rfam.pdb.gz")


    """
    mysql `python $params.rfamprod/scripts/view/mysql_options.py $params.db` < sql/ftp_rfam_pdb.sql > Rfam.pdb
    gzip Rfam.pdb
    """

}

process generate_clanin_file {
    input:
    path("Rfam.pdb.gz")
    
    output:
    val('done')


    """
    python $params.rfamprod/scripts/release/clanin_file_generator.py $params.release_ftp
    """

}


workflow generate_ftp_files {
    tree_files \
    | rfam2go_file \
    | add_header \
    | checksum \
    | validate_go \
    | generate_full_region_file \
    | generate_pdb_file \
    | generate_clanin_file
}

workflow {
    stage_rfam_live()
}
