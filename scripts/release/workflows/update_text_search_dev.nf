nextflow.enable.dsl=2
params.text_search = "$params.release/text_search"
params.xml_dumper = "$params.rfamprod/scripts/export/rfam_xml_dumper.py"
params.validate = "$params.rfamprod/scripts/validation/xml_validator.py"
params.empty = "$params.rfamprod/scripts/release/check_empty.sh"


process xml_dump {  
    memory '10GB'

    input:
    val(_flag)
    
    output:
    val('xml_done')

    """
    cd ${params.rfamprod}
    source django_settings.sh
    mkdir -p $params.text_search
    mkdir $params.text_search/families
    mkdir $params.text_search/clans
    mkdir $params.text_search/motifs
    mkdir $params.text_search/genomes
    mkdir $params.text_search/full_region   
    python $params.xml_dumper --type F --out $params.text_search/families
    python $params.xml_dumper --type C --out $params.text_search/clans
    python $params.xml_dumper --type M --out $params.text_search/motifs
    python $params.xml_dumper --type G --out $params.text_search/genomes
    python $params.xml_dumper --type R --out $params.text_search/full_region 
    """
}

process xml_validate {
    input:
    val('xml_done')

    output:
    val('validate_done')

    """ 
    python $params.validate --input $params.text_search/families --log
    python $params.validate --input $params.text_search/clans --log
    python $params.validate --input $params.text_search/motifs --log
    python $params.validate --input $params.text_search/genomes --log
    python $params.validate --input $params.text_search/full_region --log
    """
}

process check_error_logs_are_empty { 
    input:
    val('validate_done')

    output:
    val('done')

    """ 
    bash $params.empty "$params.text_search/families/error.log"
    bash $params.empty "$params.text_search/clans/error.log"
    bash $params.empty "$params.text_search/motifs/error.log"
    bash $params.empty "$params.text_search/genomes/error.log"
    bash $params.empty "$params.text_search/full_region/error.log"
    """
}
process create_release_note {
    publishDir "$params.text_search", mode: "copy"
    
    input:
    val('done')
    
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
    cp $query /nfs/production/xfam/rfam/search_dumps/rfam_dev
    unlink /nfs/production/xfam/rfam/search_dumps/rfam_dev
    ln -s $params.text_search /nfs/production/xfam/rfam/search_dumps/rfam_dev
    """
}

workflow text_search {
    take: start
    emit: done
    main:
        start | xml_dump \
        | xml_validate \
        | check_error_logs_are_empty \
        | create_release_note \
        | index_data_on_dev \
        | set { done }
    
    
}
workflow {
    text_search(Channel.of('start'))
}