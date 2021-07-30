nextflow.enable.dsl=2

process generate_clan_files {

    output:
    path('*.txt')

    """
    mkdir ~/releaseX/clan_competition
    mkdir ~/releaseX/clan_competition/sorted  
    python $baseDir/scripts/release/clan_file_generator.py --dest-dir $baseDir/releaseX/clan_competition --clan-acc CL00001 --cc-type PDB
    """

}

process sort_clan_files {
    input:
    path(query)

    output:
    path('*.txt')

    """
    cd $baseDir/releaseX/clan_competition
    for file in ./CL*; do sort -k2 -t $'\t' ${file:2:7}.txt > sorted/${file:2:7}_s.txt; done
    """

}

process run_clan_competition { 
    
    input:
    path(query)

    output:
    path('clan_file_sorted.txt')

    """
    python scripts/processing/clan_competition.py --input releaseX/clan_competition/sorted --pdb
    """
}

process get_new_families {

    """
    python pdb_new_fam_db.py
    """
}



workflow clan_compete {

    generate_clan_files \
    | sort_clan_files \
    | run_clan_competition

}