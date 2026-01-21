process apicuron_report {

    input:
    val(_flag)

    output:
    path('bulk_report_svn.json')

    """
    python $params.rfamprod/scripts/apicuron/extract_svn_info.py --end-rev $endrev --start-rev $startrev --svn $svn
    python $params.rfamprod/scripts/apicuron/bulk_upload.py -f bulk_report_svn.json
    """
}

workflow apicuron {
  take: start
  emit: done
  main:
    start
    | apicuron_report \
    | set { done }
}

workflow {
  apicuron(Channel.of('start'))
}
