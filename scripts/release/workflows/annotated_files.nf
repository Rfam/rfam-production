nextflow.enable.dsl=2

params.rfamprod = "/nfs/production/xfam/users/rfamprod/code/rfam-production"
params.release = "/hps/nobackup/production/xfam/rfam/RELEASES/14.7/"
params.release_ftp = "/hps/nobackup/production/xfam/rfam/RELEASES/14.7/ftp"

process generate_seed_files {
    memory '10GB'
    
    output:
    path("Rfam.seed")

    """
    rm -rf ${params.release_ftp}/seed/
    mkdir -p ${params.release_ftp}/seed/
    python ${params.rfamprod}/scripts/export/generate_ftp_files.py --acc all --seed --dest-dir ${params.release_ftp}/seed/
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
    perl seed-desc-to-cm.pl Rfam.seed Rfam.nodesc.cm > Rfam.cm
    """

}
process generate_archive_zip { 
    publishDir "${params.release_ftp}/cm", mode: "copy"

    input:
    path(query)

    output:
    path("Rfam.cm.gz")

    """
    cmstat $query > cmstat_file.txt
    python scripts/release/rfam_cm_check.py --cm-file $query --stat-file cmstat_file.txt
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
    rm -f RF0*.cm
    grep ACC $query | sed -e 's/ACC\\s\\+//g' | sort | uniq > list.txt
    cmfetch --index $query
    while read p; do echo ${p}; cmfetch $query ${p} > "${p.cm}" ; done <list.txt
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
    python load_cm_seed_in_db.py ${params.release_ftp}/seed/Rfam.seed ${params.release_ftp}/cm/Rfam.cm
    """

}

workflow generate_annotated_files {
    emit: 
        rfam_cm
    main:
        generate_seed_files \
        | generate_cm_files \
        | rewrite_cm_file | set { rfam_cm }
        rfam_cm \
        | generate_archive_zip
        rfam_cm \
        | create_tar_file \
        | load_into_rfam_live
}

workflow {
    generate_annotated_files()
}
