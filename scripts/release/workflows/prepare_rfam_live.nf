nextflow.enable.dsl=2

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
    perl ${params.perl_path}/jiffies/release/make_rfam_keywords_table.pl
    """
}
process update_taxonomy_websearch { 
    memory '10GB'
    
    input:
    val('keyword_complete')
    
    output:
    val('update_complete')
    
    """
    perl ${params.perl_path}/jiffies/updateTaxonomyWebsearch.pl
    """
}

workflow prepare_rfam_live {
    emit:
        done
    populate_rfam_live \
    | make_keywords \
    | update_taxonomy_websearch
    | set { done }
    
}

workflow {
    prepare_rfam_live()
}
