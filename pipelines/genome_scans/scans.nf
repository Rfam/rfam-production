process scan {
    memory { ["RF00002", "RF00005", "RF00177", "RF02542"].contains(family) ? "4GB" : "10GB" }
    queue 'short'

    input:
    val(family)

    output:
    path(family)

    """
    rfco.pl $family
    cd $family
    rfsearch.pl
    rfmake.pl
    cd ..
    """
}

process check_in {
    input:
    path(family)

    """
    rfci.pl $family -m "updated with rfamseq15"
    """
}

workflow {
    Channel.fromPath(params.families)
    | splitText(by: params.chunkSize)
    | map{it -> it.trim()}
    | scan
    | check_in
}
