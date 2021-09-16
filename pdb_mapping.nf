nextflow.enable.dsl=2

process setup_files {
    publishDir "$baseDir", mode: "copy"

    output:
    path("pdb_seqres.txt")

    """
    rm -f $baseDir/PDB_RFAM_X_Y.tbl
    rm -f $baseDir/pdb_seqres.txt
    rm -rf $baseDir/Rfam.cm*
    wget http://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
    gunzip Rfam.cm.gz
    mv Rfam.cm $baseDir
    cmpress $baseDir/Rfam.cm
    wget ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz
    gunzip pdb_seqres.txt.gz
    """
}

process remove_illegal_characters {  
    input:
    path(query)

    output:
    path("pdb_trimmed_noillegals.fa")

    """ 
    sed -e '/^[^>]/s/[^ATGCURYMKSWHBVDatgcurymkswhbvd]/N/g' $query > pdb_trimmed_noillegals.fa
    """
}

process run_cmscan {
    memory '10GB'

    input:
    path(query)
    
    output:
    path('*.tbl')
    
    """
    cmscan -o ${query}.output --tblout ${query}.tbl --cut_ga $baseDir/Rfam.cm $query
    """
}

process combine_cmscan_results {
    publishDir "$baseDir", mode: "copy"
    
    input:
    path('raw*.tbl')
    
    output:
    path('PDB_RFAM_X_Y.tbl')
    
    """
    grep -v '#' raw*.tbl | sort > PDB_RFAM_X_Y.tbl 
    """
}

process create_text_file_for_db {
    publishDir "$baseDir", mode: "copy"

    input:
    path(query)
    
    output:
    path('pdb_full_region_*.txt')
    
    """
    python $baseDir/scripts/processing/infernal_2_pdb_full_region.py --tblout $query --dest-dir .
    """
}

process import_db_and_generate_clan_files {
    input:
    path(query)
    
    output:
    path('CL*.txt')

    """
    python $baseDir/pdb_mapping/pdb_full_region_table.py --file $query
    mkdir -p $baseDir/releaseX/clan_competition/sorted  
    python $baseDir/scripts/release/clan_file_generator.py --dest-dir . --cc-type PDB
    """
}

process sort_clan_files {
    publishDir "$baseDir/releaseX/clan_competition/sorted", mode: "copy"
    
    input:
    path(query)

    output:
    path('*_s.txt')

    """ 
    for file in ./CL*; do sort -k2 -t \$'\t' \${file:2:7}.txt > \${file:2:7}_s.txt; done
    """
}

process run_clan_competition { 
    publishDir "$baseDir/releaseX/clan_competition", mode: "copy"

    input:
    path(query)

    output:
    path('*')

    """
    python $baseDir/scripts/processing/clan_competition.py --input $baseDir/releaseX/clan_competition/sorted --pdb
    """
}

process get_new_families {
    input:
    path(query)
    output:
    path('pdb_families.txt')

    """
    python $baseDir/pdb_mapping/pdb_families.py > $baseDir/pdb_families.txt
    """
}

process update_ftp {
    input:
    path(query)

    output:
    val('done')

    """
    rm /nfs/ftp/pub/databases/Rfam/.preview/pdb_full_region.txt.gz
    cp $query /nfs/ftp/pub/databases/Rfam/.preview/pdb_full_region.txt
    gzip /nfs/ftp/pub/databases/Rfam/.preview/pdb_full_region.txt
    """
}


process create_validate_xml_families {
    errorStrategy 'finish'

    input:
    path(query)

    output:
    path('relX_text_search/families/error.log')

    """
    source $baseDir/django_settings.sh
    rm -rf $baseDir/relX_text_search/families
    mkdir -p $baseDir/relX_text_search/families
    python $baseDir/scripts/export/rfam_xml_dumper.py --type F --out relX_text_search/families
    python $baseDir/scripts/validation/xml_validator.py --input relX_text_search/families --log
    bash $baseDir/pdb_mapping/check_empty.sh
    """
}

process index_data_on_rfam_dev {
    input:
    path(query)

    output:
    val('done')

    """
    cd_main && cd search_dumps
    rm -rf /nfs/production/xfam/rfam/search_dumps/rfam_dev/families/
    cp -r /nfs/production/xfam/users/rfamprod/code/rfam-production/relX_text_search/families/ rfam_dev
    """

}

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

workflow pdb_mapping {
    emit:
        pdb_txt
        new_families
    main:
    setup_files \
    | splitFasta( record: [id: true, desc:true, text: true] ) \
    | filter { record -> record.desc =~ /^mol:na.*/ } \
    | collectFile( name:"pdb_trimmed_noprot.fa") {
        it.text
    } \
    | remove_illegal_characters \
    | splitFasta ( by:300, file:true ) \
    | run_cmscan \
    | collect \
    | combine_cmscan_results \
    | create_text_file_for_db | set { pdb_txt }
    pdb_txt
    | import_db_and_generate_clan_files \
    | sort_clan_files \
    | collect \
    | run_clan_competition \
    | get_new_families | set {new_families}
}

workflow ftp {
    take: pdb_txt
    main:
    pdb_txt \
    | update_ftp
}

workflow update_search_index {
    take: new_families
    main:
    new_families \
    | create_validate_xml_families \
    | index_data_on_rfam_dev
}

workflow update_website_db {
    take: pdb_txt
    main:
    pdb_txt \
    | sync_db
}

workflow {
    pdb_mapping()
    ftp(pdb_mapping.out.pdb_txt)
    update_search_index(pdb_mapping.out.new_families)
    update_website_db(pdb_mapping.out.pdb_txt)
}

workflow.onComplete = {

def process = ["python", "pdb_mapping/send_notification.py"].execute()
process.waitFor()
println process.err.text
println process.text

println "--------------------------"
println "PDB Mapping Pipeline Summary"
println "--------------------------"
println "Started at  : ${workflow.start}"
println "Completed at: ${workflow.complete}"
println "Duration    : ${workflow.duration}"
println "Success     : ${workflow.success}"
println "workDir     : ${workflow.workDir}"
println "exit status : ${workflow.exitStatus}"

}