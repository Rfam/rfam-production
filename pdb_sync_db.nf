nextflow.enable.dsl=2

include { pdb_mapping } from './pdb_mapping'

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

workflow update_website_db {
    take: pdb_txt
    main:
    pdb_txt \
    | sync_db
}

workflow {
    update_website_db(pdb_mapping.out.pdb_txt)
}
