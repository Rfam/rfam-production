nextflow.enable.dsl=2
params.text_search = "$params.release/text_search"
params.xml_dumper = "$params.rfamprod/scripts/export/rfam_xml_dumper.py"
params.validate = "$params.rfamprod/scripts/validation/xml_validator.py"
params.empty = "$params.rfamprod/scripts/release/check_empty.sh"

// Validate required parameters at startup
if (!params.rfamprod) {
    error "Missing required parameter: params.rfamprod"
}
if (!params.release) {
    error "Missing required parameter: params.release"
}
if (!params.releasex) {
    error "Missing required parameter: params.releasex"
}
if (!params.rfam_dev) {
    error "Missing required parameter: params.rfam_dev"
}


process xml_dump {  
    memory '32GB'
    cpus 1
    errorStrategy 'terminate'

    input:
    val flag

    output:
    val true

    script:
    """
    #!/bin/bash
    set -euo pipefail

    # Set UTF-8 encoding to handle non-ASCII characters
    export LC_ALL=en_US.UTF-8
    export LANG=en_US.UTF-8
    export PYTHONIOENCODING=utf-8

    cd ${params.rfamprod} && source django_settings.sh

    #mkdir -p $params.text_search
    #mkdir $params.text_search/families
    #mkdir $params.text_search/clans
    #mkdir $params.text_search/motifs
    #mkdir $params.text_search/genomes
    #mkdir $params.text_search/full_region
    #
    #python $params.xml_dumper --type F --out $params.text_search/families
    #python $params.xml_dumper --type C --out $params.text_search/clans
    #python $params.xml_dumper --type M --out $params.text_search/motifs
    #python $params.xml_dumper --type G --out $params.text_search/genomes
    #python $params.xml_dumper --type R --out $params.text_search/full_region 

    # create directory structure
    mkdir -p ${params.text_search}/{families,clans,motifs,genomes,full_region}
    
    # run XML dumps with error checking
    python ${params.xml_dumper} --type F --out ${params.text_search}/families
    echo "Families completed"
    python ${params.xml_dumper} --type C --out ${params.text_search}/clans
    echo "Clans completed"
    python ${params.xml_dumper} --type M --out ${params.text_search}/motifs
    echo "Motifs completed"
    python ${params.xml_dumper} --type G --out ${params.text_search}/genomes
    echo "Genomes completed"
    python ${params.xml_dumper} --type R --out ${params.text_search}/full_region
    echo "Full_region completed"
    
    echo "XML dump completed successfully"
    """
}

process xml_validate {
    memory '8 GB'
    cpus 2
    errorStrategy 'terminate'

    input:
    val xml_done

    output:
    val true

    script:
    """ 
    #!/bin/bash
    set -euo pipefail

    #python $params.validate --input $params.text_search/families --log
    #python $params.validate --input $params.text_search/clans --log
    #python $params.validate --input $params.text_search/motifs --log
    #python $params.validate --input $params.text_search/genomes --log
    #python $params.validate --input $params.text_search/full_region --log

    # Validate XML files with error checking
    python ${params.validate} --input ${params.text_search}/families --log || exit 1
    python ${params.validate} --input ${params.text_search}/clans --log || exit 1
    python ${params.validate} --input ${params.text_search}/motifs --log || exit 1
    python ${params.validate} --input ${params.text_search}/genomes --log || exit 1
    python ${params.validate} --input ${params.text_search}/full_region --log || exit 1

    echo "Validation completed successfully"
    """
}

process check_error_logs_are_empty { 
    input:
    val validate_done

    output:
    val true

    script:
    """ 
    #!/bin/bash
    set -euo pipefail

    # Check that all error logs are empty
    bash ${params.empty} "${params.text_search}/families/error.log" || exit 1
    bash ${params.empty} "${params.text_search}/clans/error.log" || exit 1
    bash ${params.empty} "${params.text_search}/motifs/error.log" || exit 1
    bash ${params.empty} "${params.text_search}/genomes/error.log" || exit 1
    bash ${params.empty} "${params.text_search}/full_region/error.log" || exit 1
    
    echo "All error logs are empty"
    """
}

process create_release_note {
    publishDir "${params.text_search}", mode: "copy"
    memory '2 GB'
    errorStrategy 'terminate'
    
    input:
    val check_done
    
    output:
    path "release_note.txt"

    script:
    """ 
    #!/bin/bash
    set -euo pipefail

    #python $params.rfamprod/scripts/release/create_release_note.py --version=$params.releasex --entries=`grep -r 'entry id' $params.text_search | wc -l`
    ENTRIES=\$(grep -r 'entry id' $params.text_search | wc -l || echo "0")
    python ${params.rfamprod}/scripts/release/create_release_note.py \
      --version=${params.releasex} \
      --entries=\$ENTRIES > release_note.txt

    echo "Created release note with \$ENTRIES entries"
    """
}
process index_data_on_dev {
    input:
    path release_note

    output:
    val true

    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    # Safely remove existing symlink if it exists
    if [ -L ${params.rfam_dev} ] || [ -e ${params.rfam_dev} ]; then
        rm -f ${params.rfam_dev}
    fi
    
    # Create new symlink
    ln -s ${params.text_search} ${params.rfam_dev}
    
    # Copy release note
    cp ${release_note} ${params.rfam_dev}/
    
    echo "Indexed data on dev: ${params.rfam_dev}"

    #unlink $params.rfam_dev
    #ln -s $params.text_search $params.rfam_dev
    #cp $query $params.rfam_dev
    """
}

workflow text_search {
    take: 
    start
    
main:
    xml_dump(start)
    xml_validate(xml_dump.out)
    check_error_logs_are_empty(xml_validate.out)
    create_release_note(check_error_logs_are_empty.out)
    index_data_on_dev(create_release_note.out)

emit:
    done = index_data_on_dev.out
}

workflow {
    text_search(Channel.of('start'))
}
