nextflow.enable.dsl=2

process setup_files {
    publishDir "$params.pdb_files", mode: "copy"

    input:
    val(_flag)

    output:
    path("pdb_seqres.txt")

    """
    rm -f $params.pdb_files/PDB_RFAM_X_Y.tbl
    rm -f $params.pdb_files/pdb_seqres.txt
    rm -rf $params.pdb_files/Rfam.cm*
    wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
    gunzip Rfam.cm.gz
    mv Rfam.cm $params.pdb_files
    cmpress $params.pdb_files/Rfam.cm
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
    cmscan -o ${query}.output --tblout ${query}.tbl --cut_ga $params.pdb_files/Rfam.cm $query
    """
}

process combine_cmscan_results {
    publishDir "$params.pdb_files", mode: "copy"
    
    input:
    path('raw*.tbl')
    
    output:
    path('PDB_RFAM_X_Y.tbl')
    
    """
    grep -v '#' raw*.tbl | sort > PDB_RFAM_X_Y.tbl 
    """
}

process create_text_file_for_db {
    publishDir "$params.pdb_files", mode: "copy"

    input:
    path(query)
    
    output:
    path('pdb_full_region_*.txt')
    
    """
    python $params.rfamprod/scripts/processing/infernal_2_pdb_full_region.py --tblout $query --dest-dir .
    """
}

process import_db_and_generate_clan_files {
    input:
    path(query)
    
    output:
    path('CL*.txt')

    """
    bash $params.pdb_scripts/check_not_empty.sh $query
    python $params.pdb_scripts/pdb_full_region_table.py --file $query
    mkdir -p $params.pdb_files/clan_competition/sorted
    python $params.rfamprod/scripts/release/clan_file_generator.py --dest-dir . --cc-type PDB
    """
}

process sort_clan_files {
    publishDir "$params.pdb_files/clan_competition/sorted", mode: "copy"
    
    input:
    path(query)

    output:
    path('*_s.txt')

    """ 
    for file in ./CL*; do sort -k2 -t \$'\t' \${file:2:7}.txt > \${file:2:7}_s.txt; done
    """
}

process run_clan_competition { 
    publishDir "$params.pdb_files/clan_competition", mode: "copy"

    input:
    path(query)

    output:
    path('*')

    """
    python $params.rfamprod/scripts/processing/clan_competition.py --input $params.pdb_files/clan_competition/sorted --pdb
    """
}

process get_new_families {
    publishDir "$params.pdb_files", mode: "copy"
    
    input:
    path(query)
    
    output:
    path('pdb_families_*.txt')

    """
    python $params.pdb_scripts/pdb_families.py
    """
}

process get_ftp_file {
    publishDir "$params.pdb_files", mode: "copy"

    input:
    path(query)

    output:
    path('pdb_full_region_updated_*.txt')

    """
    python $params.pdb_scripts/pdb_ftp_file.py --dest-dir .
    """
}

process update_ftp {
    queue 'datamover'

    input:
    path(query)

    output:
    val('done')

    """
    rm -f $params.preview/pdb_full_region.txt.gz
    cp $query $params.preview/pdb_full_region.txt
    gzip $params.preview/pdb_full_region.txt
    """
}


process create_validate_xml_families {
    errorStrategy 'finish'

    input:
    path(query)

    output:
    val('xml')

    """
    source $params.rfamprod/django_settings.sh
    rm -rf $params.pdb_files/text_search/families
    mkdir -p $params.pdb_files/text_search/families
    python $params.rfamprod/scripts/export/rfam_xml_dumper.py --type F --out $params.pdb_files/text_search/families --db rfam-rel
    python $params.rfamprod/scripts/validation/xml_validator.py --input $params.pdb_files/text_search/families --log
    bash $params.pdb_scripts/check_empty.sh "$params.pdb_files/text_search/families/error.log"
    """
}

process index_data_on_rfam_dev {
    input:
    val('xml')

    output:
    val('dev_done')

    """
    rm -rf $params.search_dumps/rfam_dev/families/
    cp -r $params.pdb_files/text_search/families/ $params.search_dumps/rfam_dev/
    """

}

process index_data_on_prod {
    input:
    val('dev_done')

    output:
    val('done')

    """
    rm -rf $params.search_dumps/current_release/families/
    cp -r $params.pdb_files/text_search/families/ $params.search_dumps/current_release/
    """

}

process sync_rel_db {
    input:
    path(query)

    output:
    val('done')

    """
    python $params.pdb_scripts/pdb_full_region_table.py --file $query --database rfam-rel -nqc
    """
}

process sync_web_production_db {
    input:
    path(query)

    output:
    val('synced')
    """
    python $params.pdb_scripts/pdb_full_region_table.py --file $query --database pg -nqc
    python $params.pdb_scripts/pdb_full_region_table.py --file $query --database fb1 -nqc
    """
}

process clan_compete_rel_web {
    input:
    val('synced')

    output:
    val('done')

    """
    python $params.rfamprod/scripts/processing/clan_competition.py --input $params.pdb_files/clan_competition/sorted --pdb --sync
    """

}

process add_all_3d {
    container 'docker://rfam/rfam-3d-seed-alignments:latest'
    errorStrategy 'finish'

    input:
    path(query)

    output:
    val('3d_done')

    """
    git clone $params.seed_repo
    cd rfam-3d-seed-alignments
    python add_3d.py all --nocache
    git add data/output/RF*
    git add pdb*
    git commit -m'Auto-commit - New data from PDB pipeline'
    git push
    git log --name-status -1 > $params.pdb_files/git_output.txt
    """
}

process update_3d_message{
    input:
    val('3d_done')

    output:
    val('done')

    """
    python $params.pdb_scripts/add_3d_output.py
    """
}


workflow pdb_mapping {
    take: start
    emit:
        pdb_txt
        new_families
    main:
        start | setup_files \
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
        pdb_txt \
        | import_db_and_generate_clan_files \
        | sort_clan_files \
        | collect \
        | run_clan_competition \
        | get_new_families | set {new_families}
}

workflow ftp {
    take: new_families
    main:
        new_families \
        | get_ftp_file \
        | update_ftp
}

workflow update_search_index {
    take: new_families
    main:
        new_families \
        | create_validate_xml_families \
        | index_data_on_rfam_dev \
        | index_data_on_prod
}

workflow sync_rel_web {
    take: 
        pdb_txt
    emit: synced 
    main:
        pdb_txt \
        | sync_rel_db
        pdb_txt \
        | sync_web_production_db | set { synced } 
}

workflow add_3d {
    take:
        pdb_txt
    emit: done
    main:
        pdb_txt \
        | add_all_3d
        | update_3d_message | set { done }
}

workflow mapping_and_updates {
    take: start
    emit: done
    main:
        pdb_mapping(start)
        ftp(pdb_mapping.out.new_families)
        update_search_index(pdb_mapping.out.new_families)
        sync_rel_web(pdb_mapping.out.pdb_txt)
        clan_compete_rel_web(sync_rel_web.out.synced)
        add_3d(pdb_mapping.out.new_families) \
        | set { done }
}

workflow {
    mapping_and_updates(Channel.of('start'))
}

workflow.onComplete = {

def process = ["python", "$params.pdb_scripts/send_notification.py"].execute()
process.waitFor()
println process.err.text
println process.text

def msg = """\
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """
        .stripIndent()

sendMail(to: $params.email, subject: 'PDB pipeline execution', body: msg)
println msg
}
