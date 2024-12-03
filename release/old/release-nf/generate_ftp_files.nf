nextflow.enable.dsl=2

workflow generate_ftp_files {
    take: start
    emit: done
    main:
        fetch_accessions \
        | splitText \
        | map { it.strip() } \
        | tree_file \
        | collect \
        | tree_tarball

        tree_files(start)
        generate_full_region_file \
        | generate_pdb_file \
        | generate_clanin_file \
        | generate_3d_seed_file \
        | ena_rfam_mapping
        | set { done }
}

workflow {
    generate_ftp_files(Channel.of('start'))
}
