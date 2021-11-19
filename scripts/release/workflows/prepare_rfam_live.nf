nextflow.enable.dsl=2

params.rfamprod = "/nfs/production/xfam/users/rfamprod/code/rfam-production"
params.release = "/hps/nobackup/production/xfam/rfam/RELEASES/14.7/"

process populate_rfam_live {
    
    output:
    path('')

    """
    python ${params.release}/scripts/release/populate_rfamlive_for_release.py --all
    """
}
process make_keywords {
    publishDir "${params.release}/clan_competition/sorted", mode: "copy"
    
    input:
    path(query)
    
    output:
    path('')

    """ 
    perl make_rfam_keywords_table.pl
    """
}
process update_taxonomy_websearch { 
    publishDir "${params.release}/clan_competition/sorted", mode: "copy"
    
    input:
    path(query)
    
    output:
    path('')
    
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
