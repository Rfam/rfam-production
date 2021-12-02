nextflow.enable.dsl=2

process generate_seed_files {
    memory '10GB'
    publishDir "${params.release_ftp}/seed", mode: "copy"
    
    output:
    path("Rfam.seed")

    """
    rm -rf ${params.release_ftp}/seed/
    mkdir -p ${params.release_ftp}/seed/
    python ${params.rfamprod}/scripts/export/generate_ftp_files.py --acc all --seed --dest-dir .
    """

}
process generate_cm_files {  
    memory '10GB'
    publishDir "${params.release_ftp}/cm", mode: "copy"
    
    input:
    path(query)

    output:
    path("Rfam.cm")

    """
    rm -rf ${params.release_ftp}/cm
    mkdir -p ${params.release_ftp}/cm/
    python ${params.rfamprod}/scripts/export/generate_ftp_files.py --acc all --cm --dest-dir .
    """
}
process rewrite_cm_file { 
    publishDir "${params.release_ftp}/cm", mode: "copy"
    
    input:
    path(query)

    output:
    path("Rfam.cm")

    """
    grep -v DESC $query > Rfam.nodesc.cm
    perl ${params.perl_path}/jiffies/seed-desc-to-cm.pl ${params.release_ftp}/seed/Rfam.seed Rfam.nodesc.cm > Rfam.cm
    """

}
process checks { 
    publishDir "${params.release_ftp}/cm", mode: "copy"

    input:
    path(query)

    output:
    path("Rfam.cm")

    """
    cmstat $query > cmstat_file.txt
    python ${params.rfamprod}/scripts/release/rfam_cm_check.py --cm-file $query --stat-file cmstat_file.txt
    """

}

process generate_archive_zip {
    publishDir "${params.release_ftp}/cm", mode: "copy"

    input:
    path(query)

    output:
    path("Rfam.cm.gz")

    """
    gzip -c $query > Rfam.cm.gz
    """

}

process create_tar_file {
    publishDir "${params.release_ftp}/cm", mode: "copy"

    input:
    path(query)

    output:
    path("Rfam.tar.gz")

    """
    cd ${params.release_ftp}/cm
    rm -f RF0*.cm
    grep ACC $query | sed -e 's/ACC\\s\\+//g' | sort | uniq > list.txt
    cmfetch --index $query
    while read p; do echo \$p; cmfetch Rfam.cm \$p > "\$p.cm" ; done <list.txt
    tar -czvf Rfam.tar.gz RF0*.cm
    rm RF0*.cm
    """
}

process load_into_rfam_live {
    input:
    path(query)

    output:
    val('done')

    """
    python ${params.rfamprod}/scripts/release/load_cm_seed_in_db.py ${params.release_ftp}/seed/Rfam.seed $query
    """

}

workflow generate_annotated_files {
    emit: 
        rfam_cm
    main:
        generate_seed_files \
        | generate_cm_files \
        | rewrite_cm_file \
        | checks | set { rfam_cm }
        rfam_cm
        | generate_archive_zip
        rfam_cm \
        | create_tar_file 
        rfam_cm \
        | load_into_rfam_live
}

workflow {
    generate_annotated_files()
}
