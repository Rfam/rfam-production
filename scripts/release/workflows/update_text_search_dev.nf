nextflow.enable.dsl=2

params.text_search = "$params.release/text_search"
params.xml_dumper = "$params.rfamprod/scripts/export/rfam_xml_dumper.py"
params.validate = "$params.rfamprod/scripts/validation/xml_validator.py"
params.empty = "$params.rfamprod/scripts/release/check_empty.sh"

groups = Channel.from('families', 'clans', 'motifs', 'genomes', 'full_region')

process xml_dump {   
    input:
    tuple val(group), val(type) 
    
    output:
    val('xml_done')

    """
    cd ${params.rfamprod}
    source django_settings.sh
    mkdir -p $params.text_search
    mkdir $params.text_search/$group  
    python $params.xml_dumper --type $type --out $params.text_search/$group
    """
}

process xml_validate {
    input:
    tuple val('xml_done'), val(group)
    
    output:
    val('validate_done')

    """ 
    python $params.validate --input $params.text_search/$group --log
    """
}

process check_error_logs_are_empty { 
    input:
    tuple val('validate_done'), val(group)
    
    output:
    val('empty')

    """ 
    bash $params.empty "$params.text_search/$group/error.log"
    """
}
process create_release_note {
    publishDir "$params.text_search", mode: "copy"
    
    input:
    val('empty')
    
    output:
    path("release_note.txt")

    """ 
    python $params.rfamprod/scripts/release/create_release_note.txt --version=$params.releasex --entries=`grep -r 'entry id' $params.text_search | wc -l`
    """
}

process index_data_on_dev {
    input:
    path(query)

    output:
    val('dev_done')

    """
    unlink /nfs/production/xfam/rfam/search_dumps/rfam_dev
    ln -s $params.text_search /nfs/production/xfam/rfam/search_dumps/rfam_dev
    """
}

workflow text_search {
    emit: 
        done
        validated
    main: 
        Channel.from( ['families', 'F'], ['clans', 'C'], ['motifs', 'M'], ['genomes', 'G'], ['full_region', 'R'] ) \
        | xml_dump | set {done}
        done, groups \
        | xml_validate | set {validated}
        validated, groups \
        | check_error_logs_are_empty \
        | create_release_note \
        | index_data_on_dev 
}

workflow {
    text_search()
}
