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
    path('rfam2go')

    """
    python $params.rfamprod/scripts/release/add_header_rfam2go.py --input $query
    """
}

process checksum {
    publishDir "$params.release_ftp", mode: "copy"
    
    input:
    path(query)
    
    output:
    path('md5.txt')

    """
    md5sum $query > md5.txt
    """
}

process validate_go {
    publishDir "$params.release_ftp", mode: "copy"
    input:
    path(query)
    
    output:
    val('go_validated')

    """
    perl $params.perl_path/qc/validate_rfam2go.pl goselect selectgo $query
    """

}

workflow generate_rfam2go {
    emit: 
    rfam2go_file
    
    main:
    rfam2go_file \
    | add_header | set {rfam2go_file}
    rfam2go_file \
    | checksum
    rfam2go_file \
    | validate_go
}

workflow {
    generate_rfam2go()
}