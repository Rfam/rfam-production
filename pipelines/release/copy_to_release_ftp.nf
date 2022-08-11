params.ftp_dir = "$params.ftp_site/$params.releasex"

process copy_to_release_ftp {
    queue 'datamover'

    output:
    val('done)')

    """
    mkdir $params.ftp_dir
    cd $params.release_ftp
    cp -r fasta_files/ $params.ftp_dir
    cp -r database_files/ $params.ftp_dir
    cp -r rfam2go/ $params.ftp_dir
    cp Rfam.full_region.gz $params.ftp_dir
    cp Rfam.clanin $params.ftp_dir
    cp Rfam.pdb.gz $params.ftp_dir
    cp tree/Rfam.seed_tree.tar.gz $params.ftp_dir
    cp cm/Rfam.cm.gz $params.ftp_dir
    cp cm/Rfam.tar.gz $params.ftp_dir
    cp seed/Rfam.seed.gz $params.ftp_dir
    cp seed/Rfam.seed_tree.tar.gz $params.ftp_dir

    cp -r $params.releases_dir/genome_browser_hub/ $params.ftp_dir
    cp $params.releases_dir/COPYING $params.ftp_dir
    cp $params.releases_dir/USERMAN $params.ftp_dir
    cp $params.releases_dir/README $params.ftp_dir
    """

}

workflow {
    copy_to_release_ftp()
}
