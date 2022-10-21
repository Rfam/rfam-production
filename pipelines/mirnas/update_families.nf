process checkout_copy_seed {  
    input:
    path(mirnas)
    
    output:
    val('copied')

    """
    python $params.mirnas/copy_seed.py --csv-input $query
    """
}

process run_searches {
    input:
    tuple path(mirnas), val(copied)
    
    output:
    val('searches_done')

    """ 
    python $params.mirnas/auto_rfsearch.py --csv-input $query
    """
}

process run_makes {
    input:
    path(mirnas), val(searches_done)
    
    output:
    val('make_done')

    """ 
    python $params.mirnas/auto_rfmake.py --csv-input $query
    """
}

process add_refs_to_desc {
    input:
    path(mirnas), val(make_done)
    
    output:
    val('refs_done')

    """ 
    python $params.mirnas/auto_addref.py --csv-input $query
    """
}

process update_desc {
    input:
    path(mirnas), val(refs_done)
    
    output:
    val('desc_done')

    """ 
    python $params.mirnas/update_desc.py --csv-input $query
    """
}

process qc_checks {
    input:
    path(mirnas), val(desc_done)
    
    output:
    val('qc_done')

    """ 
    python $params.mirnas/auto_rqc.py --csv-input $query
    """
}

workflow update_families {
    take: mirnas_csv
    main: 
        checkout_copy_seed(mirnas_csv)
        run_searches(mirnas_csv, checkout_copy_seed.out.copied)
        run_makes(mirna_csv, run_searches.out.searches_done)
        add_refs_to_desc(mirna_csv, run_makes.out.make_done)
        update_desc(mirna_csv, add_refs_to_desc.out.refs_done)
        qc_checks(mirna_csv, update_desc.out.desc_done)
    
}

workflow {
    update_families(Channel.of(params.in))
}
