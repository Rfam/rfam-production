nextflow.enable.dsl=2

process generate_seed_files {
    memory '10GB'

    output:
    val('done)')

    """
    rm -rf ${params.release_ftp}/seed/
    mkdir -p ${params.release_ftp}/seed/
    python ${params.rfamprod}/scripts/export/generate_ftp_files.py --acc all --seed --dest-dir "${params.release_ftp}/seed"
    """

}

workflow copy_to_release_ftp {
    emit:
        done
    main:
        copy_seed_files \
        copy_cm_files \
        copy_genome_browser \
        copy_fasta \
        copy_tree \
        copy_files \

        | set { done }
}

workflow {
    copy_to_release_ftp()
}
