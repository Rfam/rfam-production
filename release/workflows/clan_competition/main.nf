nextflow.enable.dsl=2

process generate_clan_files {
    memory '10GB'
    
    input:
    val(_flag)
    
    output:
    path('CL*.txt')

    script:
    """
    mkdir -p ${params.release}/clan_competition/sorted  
    python ${params.rfamprod}/scripts/release/clan_file_generator.py --dest-dir . --cc-type FULL
    """
}

process sort_clan_files {
    publishDir "${params.release}/clan_competition/sorted", mode: "copy"
    
    input:
    path(query)
    
    output:
    path('*_s.txt')

    script:
    """ 
    #for file in ./CL*; do sort -k2 -t \$'\t' \${file:2:7}.txt > \${file:2:7}_s.txt; done
    # recommended fix:
    for file in CL*.txt; do 
        base=\$(basename \$file .txt)
        sort -k2 -t \$'\t' \$file > \${base}_s.txt
    done
    """
}

process run_clan_competition { 
    publishDir "${params.release}/clan_competition/sorted", mode: "copy"
    
    input:
    path(query)
    
    output:
    val('done')

    script:
    """
    python ${params.rfamprod}/scripts/processing/clan_competition.py --input ${params.release}/clan_competition/sorted --full
    """
}

workflow clan_competition {
    take: start
    main:
        done = start | generate_clan_files \
               | sort_clan_files \
               | run_clan_competition
    emit: done
}

