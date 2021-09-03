nextflow.enable.dsl=2

process update_website_db {
     
    input:
    path(query)

    output:
    path(query) 

    """
    python update_pdb_rel_db.py --file $query
    become xfm_adm
    ssh ves-oy-b7
    sudo /etc/init.d/httpd stop
    exit
    exit
    become mysql-rel-4442
    sync-mysql-fb --dc=FB1
    sudo /etc/init.d/httpd start
    """
}


workflow pdb_updates {

    update_website_db

}

workflow {
    pdb_updates()
}
