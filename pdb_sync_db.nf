nextflow.enable.dsl=2

pdb_text = file(params.in)

process sync_db {
    input:
    path(query)

    output:
    val('done')

    """
    cd_code && cd rfam-production
    python $baseDir/pdb_mapping/pdb_full_region_table.py --file $query --database rfam-rel -nqc
    become mysql-rel-4442
    yes | sync-mysql-fb --dc=FB1
    exit
    """
}

workflow {
    sync_db(pdb_text)
}
