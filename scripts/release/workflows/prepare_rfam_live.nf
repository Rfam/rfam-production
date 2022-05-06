nextflow.enable.dsl=2

process populate_rfam_live {
    input:
    val(_flag)
    
    output:
    val('populate_complete')

    """
    python $params.rfamprod/scripts/release/populate_rfamlive_for_release.py --all
    """
}
process make_keywords {
    
    input:
    val('populate_complete')
    
    output:
    val('keyword_complete')

    """ 
    perl $params.perl_path/jiffies/release/make_rfam_keywords_table.pl
    """
}
process update_taxonomy_websearch { 
    memory '10GB'
    
    input:
    val('keyword_complete')
    
    output:
    val('taxonomy_complete')
    
    """
    perl $params.perl_path/jiffies/updateTaxonomyWebsearch.pl
    """
}

process update_version_table {
    input:
    val('taxonomy_complete')
    
    output:
    val('version_complete')
    
    """
    python $params.rfamprod/scripts/release/update_version_table.py --version $params.releasex
    """

}

workflow prepare_rfam_live {
    take: start 
    emit: done
    start | populate_rfam_live \
    | make_keywords \
    | update_taxonomy_websearch \
    | update_version_table \
    | set { done }
    
}

workflow {
    prepare_rfam_live(Channel.of('start'))
}
