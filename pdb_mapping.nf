nextflow.enable.dsl=2

process setup_files {

    publishDir "$baseDir", mode: "copy"

    output:
    path("sample_seqs.txt")

    """
    rm -f $baseDir/PDB_RFAM_X_Y.tbl
    rm -f $baseDir/pdb_seqres.txt
    rm -rf $baseDir/cmscan_results   
    mkdir $baseDir/cmscan_results
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
    path(query)

    output:
    path('PDB_RFAM_X_Y.tbl')

    """
    cat *.tbl | sort | grep -v "#" > PDB_RFAM_X_Y.tbl
    """

}

workflow pdb_mapping {

    // setup_files \
    | splitFasta( record: [id: true, desc:true, text: true] ) \
    | filter { record -> record.desc =~ /^mol:na.*/ } \
    | collectFile( name:"pdb_trimmed_noprot.fa") {
        it.text
    } \
    | remove_illegal_characters \
    | splitFasta ( by:10, file:true ) \
    | run_cmscan \
    | collect \
    | combine_cmscan_results 
}


workflow {
    pdb_mapping()
}