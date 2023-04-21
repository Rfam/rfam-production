process scan {
    memory { ["RF00002", "RF00005", "RF00177", "RF02542"].contains(family) ? "4GB" : "10GB" }
    queue 'short'

    input:
    tuple val(family), val(threshold)

    output:
    path('${family}')

    """
    rfco.pl $family
    cd $family
    rfsearch.pl -dbchoice testrfamseq
    rfmake.pl
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
    Channel.fromPath(params.input)
    | splitText(by: params.chunkSize)
    | scan
    | check_in
}
