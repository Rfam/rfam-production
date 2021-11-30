nextflow.enable.dsl=2

params.rfamprod = "/nfs/production/xfam/users/rfamprod/code/rfam-production"
params.release = "/hps/nobackup/production/xfam/rfam/RELEASES/14.7"

process generate_clan_files {
    memory '10GB'
    
    output:
    path('CL*.txt')

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

    """ 
    for file in ./CL*; do sort -k2 -t \$'\t' \${file:2:7}.txt > \${file:2:7}_s.txt; done
    """
}

process run_clan_competition { 
    publishDir "${params.release}/clan_competition/sorted", mode: "copy"
    
    input:
    path(query)
    
    output:
    val('done')

    """
    python ${params.rfamprod}/scripts/processing/clan_competition.py --input ${params.release}/clan_competition/sorted --full
    """
}

workflow clan_competition {
    generate_clan_files \
    | sort_clan_files \
    | run_clan_competition
    
}

workflow {
    clan_competition()
}
