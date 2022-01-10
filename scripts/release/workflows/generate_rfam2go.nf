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
    python $params.rfamprod/scripts/release/add_header_rfam2go.py --input $query --version $params.releasex
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

workflow generate_rfam2go {
    emit:
        done
    main:
        rfam2go_file \
        | add_header \
        | checksum \
        | set { done }
}

workflow {
    generate_rfam2go()
}