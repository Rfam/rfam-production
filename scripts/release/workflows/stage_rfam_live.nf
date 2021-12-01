nextflow.enable.dsl=2

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
    python ${params.rfamprod}/scripts/release/restore_mysqldump.py --rel
    python ${params.rfamprod}/scripts/release/restore_mysqldump.py --public
    """
}


workflow stage_rfam_live {
    mysql_dump \
    | restore_mysql 
    
}

workflow {
    stage_rfam_live()
}
