nextflow.enable.dsl=2

params.rfamprod = "/nfs/production/xfam/users/rfamprod/code/rfam-production"

process populate_rfam_live {
    
    output:
    val('populate_complete')

    """
    python ${params.rfamprod}/scripts/release/populate_rfamlive_for_release.py --all
    """
}
process make_keywords {
    
    input:
    val('populate_complete')
    
    output:
    val('keyword_complete')

    """ 
    perl make_rfam_keywords_table.pl
    """
}
process update_taxonomy_websearch { 
    
    input:
    val('keyword_complete')
    
    output:
    val('update_complete')
    
    """
    perl updateTaxonomyWebsearch.pl
    """
}

workflow prepare_rfam_live {
    populate_rfam_live \
    | make_keywords \
    | update_taxonomy_websearch
    
}

workflow {
    prepare_rfam_live()
}
