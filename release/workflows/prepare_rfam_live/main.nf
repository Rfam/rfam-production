nextflow.enable.dsl=2

process populate_rfam_live {
    input:
    val _flag
    
    output:
    val 'populate_complete'

    script:
    """
    python $params.rfamprod/scripts/release/populate_rfamlive_for_release.py --all
    """
}

// the db this points to doesn't exist anymore
process make_keywords {
    
    input:
    val populate_complete
    
    output:
    val 'keyword_complete'

    script:
    """ 
    perl $params.perl_path/jiffies/release/make_rfam_keywords_table.pl
    """
}
process update_taxonomy_websearch { 
    memory '12GB'
    
    input:
    val keyword_complete
    
    output:
    val 'taxonomy_complete'
    
    script:
    """
    perl $params.perl_path/jiffies/updateTaxonomyWebsearch.pl
    """
}

process update_version_table {
    input:
    val taxonomy_complete
    
    output:
    val 'version_complete'
    
    script:
    """
    python $params.rfamprod/scripts/release/update_version_table.py --version $params.releasex
    """

}

process update_rnacentral_descriptions {
    input:
    val version_complete

    output:
    val 'rnacentral_complete'

    script:
    """
    python $params.rfamprod/scripts/release/update_rnacentral_description.py
    """
}

workflow prepare_rfam_live {
    take: 
        start 
    
    main:
        pop_result = populate_rfam_live(start)
        key_result = make_keywords(pop_result)
        tax_result = update_taxonomy_websearch(key_result)
        ver_result = update_version_table(tax_result)
        done = update_rnacentral_descriptions(ver_result)
    
    emit:
        done
}

workflow {
    prepare_rfam_live(Channel.of('start'))
}
