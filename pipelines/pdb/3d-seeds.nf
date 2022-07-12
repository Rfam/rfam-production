process add_all_3d {
    container 'docker://rfam/rfam-3d-seed-alignments:latest'

    output:
    val('done')

    """
    git clone $params.3d_repo
    add_3d.py all --nocache
    """

}
workflow {
    add_all_3d()
}
