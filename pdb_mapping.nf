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
    $baseDir/cmpress $baseDir/Rfam.cm
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
    $baseDir/cmscan -o ${query}.output --tblout ${query}.tbl --cut_ga $baseDir/Rfam.cm $query
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
    python $baseDir/utils/pdb_full_region_table.py --file $query
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

    publishDir "$baseDir", mode: "copy"

    input:
    path(query)

    """
    python $baseDir/scripts/release/pdb_families.py > pdb_families.txt
    """
}

workflow pdb_mapping {
    
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
    | create_text_file_for_db \
    | import_db_and_generate_clan_files \
    | sort_clan_files \
    | collect \
    | run_clan_competition \
    | get_new_families

}

workflow {
    pdb_mapping()
}