nextflow.enable.dsl=2

process mysql_dump {
    input:
    val(_flag)

    output:
    val('mysqldump_done')

    """
    mysqldump `python $params.rfamprod/scripts/view/mysql_options.py $params.db` --single-transaction --add-locks --lock-tables --add-drop-table --dump-date --comments --allow-keywords --max-allowed-packet=1G --no-create-db rfam_live > $params.release/rfam_live_rel_${params.releasex}.sql
    sed -i.bak '/USE `rfam_live`;/d' $params.release/rfam_live_rel_${params.releasex}.sql
    """
}

process restore_mysql {

    input:
    val('mysqldump_done')

    output:
    val('done')

    """
    mysql `python $params.rfamprod/scripts/view/mysql_options.py $params.db_rel` <<< "Create schema $params.db_schema_name;Use $params.db_schema_name;source $params.release/rfam_live_rel_${params.releasex}.sql;"
    mysql `python $params.rfamprod/scripts/view/mysql_options.py $params.db_pub` <<< "Create schema $params.db_schema_name;Use $params.db_schema_name;source $params.release/rfam_live_rel_${params.releasex}.sql;"
    """
}


workflow stage_rfam_live {
    take: start
    emit: done
    main:
        start | mysql_dump \
        | restore_mysql \
        | set { done }

}

workflow {
    stage_rfam_live(Channel.of('start'))
}
