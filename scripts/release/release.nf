nextflow.enable.dsl=2

process generate_seed_files {
    output:
    path("Rfam.seed")

    """
    python $baseDir/generate_ftp_files.py --acc all --seed --dest-dir /hps/nobackup/production/xfam/rfam/RELEASES/14.7/release/ftp/seed/
    """

}
process generate_cm_files {  
    input:
    path(query)

    output:
    path("Rfam.cm")

    """
    python $baseDir/generate_ftp_files.py --acc all --cm --dest-dir /hps/nobackup/production/xfam/rfam/RELEASES/14.7/release/ftp/cm/
    """
}
process rewrite_cm_file { 
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
    publishDir "$baseDir", mode: "copy"

    input:
    path(query)

    output:
    path("Rfam.cm.gz")

    // convert this to a check 
    """
    # check that Rfam.cm contains the correct number of families
    cmstat $query | grep -v '#' | wc -l

    # check the number of DESC lines - should be 2 * number of families
    grep DESC $query | wc -l

    gzip -c $query > Rfam.cm.gz
    """

}
process create_tar_file {
    publishDir "$baseDir", mode: "copy"

    input:
    path(query)

    output:
    path("Rfam.tar.gz")

    """
    rm -f RF0*.cm
    grep ACC $query | sed -e 's/ACC\s\+//g' | sort | uniq > list.txt
    cmfetch --index $query
    while read p; do echo ${p}; cmfetch $query ${p} > "${p.cm}" ; done <list.txt
    tar -czvf Rfam.tar.gz RF0*.cm
    rm RF0*.cm
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
        | create_tar_file
}
