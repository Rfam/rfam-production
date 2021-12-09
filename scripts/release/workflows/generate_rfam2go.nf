nextflow.enable.dsl=2

process rfam2go_file {
    output:
    path('rfam2go')

    """
    perl $params.perl_path/export/rfam2go.pl
    """
}

process add_header {
    publishDir "$params.release_ftp", mode: "copy"
    
    input:
    path(query)
    
    output:
    val('rfam2go_header')

    """
    python $params.rfamprod/scripts/release/add_header_rfam2go.py --input $params.release_ftp/rfam2go
    """
}

process checksum {
    publishDir "$params.release_ftp", mode: "copy"
    
    input:
    val('rfam2go_header')
    
    output:
    path('md5.txt')

    """
    md5sum $params.release_ftp/rfam2go > md5.txt
    """
}

process validate_go {
    input:
    path(query)
    
    output:
    val('go_validated')

    """
    perl $params.perl_path/qc/validate_rfam2go.pl goselect selectgo $params.release_ftp/rfam2go
    """

}

workflow generate_rfam2go {
    rfam2go_file \
    | add_header \
    | checksum \
    | validate_go
}

workflow {
    generate_rfam2go()
}