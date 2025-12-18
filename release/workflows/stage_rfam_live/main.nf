process mysql_dump {
    errorStrategy 'retry'
    maxRetries 2
    memory { 8.GB * task.attempt }

    input:
        val(_flag)
    
    output:
        val('mysqldump_done')
    
    script:
    """
    # mysqldump `python $params.rfamprod/scripts/view/mysql_options.py $params.db` --single-transaction --add-locks --lock-tables --add-drop-table --dump-date --comments --allow-keywords --max-allowed-packet=1G --no-create-db rfam_live > $params.release/rfam_live_rel_${params.releasex}.sql
    # sed -i.bak '/USE `rfam_live`;/d' $params.release/rfam_live_rel_${params.releasex}.sql
    mysqldump \
        --host=${params.db.host} \
        --port=${params.db.port} \
        --user=${params.db.user} \
        --password=${params.db.password} \
        --single-transaction \
        --add-locks \
        --lock-tables \
        --add-drop-table \
        --dump-date \
        --comments \
        --allow-keywords \
        --max-allowed-packet=1G \
        --no-create-db \
        ${params.db.name} \
        > ${params.release}/rfam_live_rel_${params.releasex}.sql
    
    sed -i.bak '/USE `rfam_live`;/d' ${params.release}/rfam_live_rel_${params.releasex}.sql
    """
}

process restore_mysql {
    errorStrategy 'retry'
    maxRetries 2
    memory { 8.GB * task.attempt }
    time { 4.h * task.attempt }
    
    input:
        val(mysqldump_done)
    
    output:
        val('done')
    
    script:
    """
    #mysql `python $params.rfamprod/scripts/view/mysql_options.py $params.db_rel` <<< "Create schema $params.db_schema_name;Use $params.db_schema_name;source $params.release/rfam_live_rel_${params.releasex}.sql;"
    #mysql `python $params.rfamprod/scripts/view/mysql_options.py $params.db_pub` <<< "Create schema $params.db_schema_name;Use $params.db_schema_name;source $params.release/rfam_live_rel_${params.releasex}.sql;"

    # Restore to release database (db_rel)
    mysql \
        --host=${params.db_rel.host} \
        --port=${params.db_rel.port} \
        --user=${params.db_rel.user} \
        --password=${params.db_rel.password} \
        -e "CREATE DATABASE IF NOT EXISTS ${params.db_schema_name};"
    
    mysql \
        --host=${params.db_rel.host} \
        --port=${params.db_rel.port} \
        --user=${params.db_rel.user} \
        --password=${params.db_rel.password} \
        ${params.db_schema_name} < ${params.release}/rfam_live_rel_${params.releasex}.sql

    # Restore to public database (db_pub)
    mysql \
        --host=${params.db_pub.host} \
        --port=${params.db_pub.port} \
        --user=${params.db_pub.user} \
        --password=${params.db_pub.password} \
        -e "CREATE DATABASE IF NOT EXISTS ${params.db_schema_name};"
    
    mysql \
        --host=${params.db_pub.host} \
        --port=${params.db_pub.port} \
        --user=${params.db_pub.user} \
        --password=${params.db_pub.password} \
        ${params.db_schema_name} < ${params.release}/rfam_live_rel_${params.releasex}.sql
    """
}


workflow stage_rfam_live {
    take: 
        start    
    main:
        start \
        | mysql_dump \
        | restore_mysql \
        | set { done }
        
    emit: 
        done 
}

workflow {
    stage_rfam_live(Channel.of('start'))
}
