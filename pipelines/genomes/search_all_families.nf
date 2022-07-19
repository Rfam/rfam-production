
process check_out_and_copy {
    publishDir "$params.test_rfamseq", mode: "copy"

    input:
    val(family)

    output:
    path('RF*')

    """
    rfco.pl $family
    cp $family $params.test_rfamseq/copies/
    """

}
process run_searches {
    publishDir "$params.test_rfamseq", mode: "copy"

    input:
    path(family)

    output:
    path('search.err')

    """
    bsub -q short -o $family/search.out -e $family/search.err "rfsearch.pl -t 25 --scpu 0 -cnompi -ignoresm"
    """
}

process run_makes {

    input:
    tuple path(family), val(threshold)

    """
    bsub -q short -o $family/make.out -e $family/make.err "rfmake.pl -t $threshold -forcethr -a -local"
    """
}

workflow {
    Channel.fromPath(params.file)
    .splitText()
    .map { file(it)}
    .view { "value: $it" } \
    | check_out_and_copy
    | run_searches
}