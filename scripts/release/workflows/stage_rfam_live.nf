nextflow.enable.dsl=2

process mysql_dump {
    
    output:
    val('mysqldump_done')

    """
    mysqldump `python $params.rfamprod/scripts/view/mysql_options.py $params.db` --single-transaction --add-locks --lock-tables --add-drop-table --dump-date --comments --allow-keywords --max-allowed-packet=1G rfam_live > $params.release_ftp/rfam_live_rel_${params.releasex}.sql
    """
}

process restore_mysql {
    
    input:
    val('mysqldump_done')
    
    output:
    val('done')

    """ 
    mysql `python $params.rfamprod/scripts/view/mysql_options.py $params.db_rel` <<< "Create schema $params.db_schema_name;Use $params.db_schema_name;source $params.release_ftp/rfam_live_rel_${params.releasex}.sql;"
    mysql `python $params.rfamprod/scripts/view/mysql_options.py $params.db_pub` <<< "Create schema $params.db_schema_name;Use $params.db_schema_name;source $params.release_ftp/rfam_live_rel_${params.releasex}.sql;"
    """
}


workflow stage_rfam_live {
    mysql_dump \
    | restore_mysql 
    
}

workflow {
    stage_rfam_live()
}
