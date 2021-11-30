nextflow.enable.dsl=2

params.rfamprod = "/nfs/production/xfam/users/rfamprod/code/rfam-production"
params.release = "/hps/nobackup/production/xfam/rfam/RELEASES/14.7"
params.releasex = "14.7"

process mysql_dump {
    
    output:
    path("rfam_live_rel${params.releasex}.sql")

    """
    cd ${params.release}
    bsub -o mysqldump.out -e mysqldump.err "mysqldump -u admin -h mysql-rfam-live -P 4445 --single-transaction --add-locks --lock-tables --add-drop-table --dump-date --comments --allow-keywords --max-allowed-packet=1G rfam_live > rfam_live_rel${params.releasex}.sql"
    """
}

process restore_mysql {
    
    input:
    path(query)
    
    output:
    val('done')

    """ 
    cd ${params.release}/rfam_live_rel${params.releasex}.sql
    python restore_mysqldump.py --rel
    python restore_mysqldump.py --public
    """
}


workflow stage_rfam_live {
    mysql_dump \
    | restore_mysql 
    
}

workflow {
    stage_rfam_live()
}
