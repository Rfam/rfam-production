process add_all_3d {
    container ''

    output:
    val('done')

    """
    add_3d.py all --nocache
    """

}
workflow ftp {
    add_all_3d()
}
