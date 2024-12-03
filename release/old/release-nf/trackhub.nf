process map_genomes {
    input:
    val(_start)

    output:
    path('genomes.csv')

    """
    curl
    """
}

process build_trackhub {
}

workflow trackhub {
    take:
        start
    emit:
        done
    main:
        start \
        | map_genomes \
        | splitCsv \
        | map { it.trim() }
        | build_trackhub \
        | collect \
        | merge
}

workflow {
    trackhub(Channel.of('start'))
}
