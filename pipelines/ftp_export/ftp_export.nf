process dirs_fetch_families {

    output:
    file('families')

    """
    rm -rf $params.ftp_exports/seed
    mkdir -p $params.ftp_exports/seed
    rm -rf $params.ftp_exports/cm
    mkdir -p $params.ftp_exports/cm
    mysql -s `python $params.scripts/view/mysql_options.py $params.db` <<< "select rfam_acc from family" > families
    """
}

process generate_seed_files {
    publishDir "$params.ftp_exports/seed", mode: "copy"
    queue 'short'

    input:
    val(acc)

    output:
    path('seed_done')

    """
    python $params.scripts/export/generate_ftp_files.py --seed --dest-dir . --acc $acc
    """

}

process generate_cm_files {
    queue 'short'
    publishDir "$params.ftp_exports/cm", mode: "copy"

    input:
    val(acc)

    output:
    val('cm_done')

    """
    python $params.scripts/export/generate_ftp_files.py --cm --dest-dir . --acc $acc
    """
}

process combined_files {

    input:
    tuple val('_flag'), val('_flag')

    output:
    tuple path('Rfam.seed'), path('Rfam.cm')

    """
    python $params.scripts/export/generate_ftp_files.py --combine --seed --dest-dir $params.ftp_exports/seed
    python $params.scripts/export/generate_ftp_files.py --combine --cm --dest-dir $params.ftp_exports/cm
    """

}

process rewrite_cm_file {
    publishDir "$params.ftp_exports/cm", mode: "copy"

    input:
    val('cm_done')

    output:
    path("Rfam.cm")

    """
    grep -v DESC $params.ftp_exports/cm/Rfam.cm > Rfam.nodesc.cm
    perl $params.perl_path/jiffies/seed-desc-to-cm.pl $params.ftp_exports/seed/Rfam.seed Rfam.nodesc.cm > Rfam.cm
    """

}
process checks {
    publishDir "$params.ftp_exports/cm", mode: "copy"

    input:
    path(query)

    output:
    path("Rfam.cm")

    """
    cmstat $query > cmstat_file.txt
    python $params.scripts/release/rfam_cm_check.py --cm-file $query --stat-file cmstat_file.txt
    """

}

process generate_archive_zip {
    publishDir "$params.ftp_exports/cm", mode: "copy"

    input:
    path(query)

    output:
    path("Rfam.cm.gz")

    """
    gzip -c $query > Rfam.cm.gz
    """

}

process create_tar_file {
    publishDir "$params.ftp_exports/cm", mode: "copy"

    input:
    path(query)

    output:
    val("done")

    """
    cd $params.ftp_exports/cm
    mkdir cm_files
    mv *RF0* cm_files
    cmfetch --index $query
    while read p; do echo \$p; cmfetch Rfam.cm \$p > "\$p.cm" ; done <list.txt
    tar -czvf Rfam.tar.gz RF0*.cm
    rm RF0*.cm
    """
}

process tree_files {
    memory '10GB'

    input:
    path(query)

    output:
    val('tree_done')

    """
    rm -rf $params.ftp_exports/tree
    mkdir -p $params.ftp_exports/tree
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

    input:
    path(query)

    output:
    val('done')

    """
    python $params.scripts/release/clanin_file_generator.py $params.ftp_exports
    """
}

process copy_to_preview {
    queue 'datamover'

    input:
    val('done')

    output:
    val('done')

    """
    cp $params.ftp_exports/Rfam.clanin $params.ftp_preview
    cp $params.ftp_exports/cm/Rfam.cm.gz $params.ftp_preview
    cp $params.ftp_exports/Rfam.full_region.gz $params.ftp_preview
    cp $params.ftp_exports/Rfam.pdb.gz $params.ftp_preview
    cp $params.ftp_exports/seed/Rfam.seed.gz $params.ftp_preview
    """

}



workflow generate_ftp_files {
    emit:
        rfam_cm
        done

    main:
        dirs_fetch_families | set { families }
        families \
        | splitText \
        | generate_seed_files | set { seed }
        families \
        | splitText \
        | generate_cm_files | set { cm }
        seed, cm \
        | combined_files | set { seed_file, cm_file }
        | rewrite_cm_file \
        | checks | set { rfam_cm }
        cm_file \
        | generate_archive_zip
        cm_file \
        | create_tar_file
        seed_file \
        | tree_files
        generate_full_region_file \
        | generate_pdb_file \
        | generate_clanin_file
        | set { done }
        done \
        | copy_to_preview | set { done }
}

workflow {
    generate_ftp_files()
}
