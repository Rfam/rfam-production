process apicuron_report {
    queue 'datamover'

    input:
    val(ready)

    output:
    path('bulk_report_svn.json')

    script:
    """
    # where are endrev, startrev, svn defined?
    # assuming we do not need to run this anymore
    # python $params.rfamprod/scripts/apicuron/extract_svn_info.py --end-rev \$endrev --start-rev \$startrev --svn \$svn

    # python $params.rfamprod/scripts/apicuron/bulk_upload.py -f bulk_report_svn.json
    # startrev: svn log -r "{2023-08-17}:HEAD" --limit 1
    # do we need to run this in the datamover queue?
    python ${projectDir}/../scripts/apicuron/bulk_upload.py \
        --start $params.svn.startrev \
        --end $params.svn.endrev \
        --svn $params.svn.svnrepo
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
