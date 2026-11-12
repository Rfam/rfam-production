/*
APICURON submission: build curation reports from the Rfam SVN history and
upload the ones APICURON does not have yet (scripts/apicuron/bulk_upload.py).

Every bulk_upload.py option comes from task.ext.args, set in nextflow.config
(process withName: 'APICURON_REPORT'). See scripts/apicuron/docs/usage.md.

Exit codes: 0 fine; 1 failure; 2 some ORCIDs need a correction in
scripts/apicuron/conf.py (the valid reports are uploaded, but the task fails
so the warning is not missed).
*/

process APICURON_REPORT {
    queue 'datamover'

    input:
    val(ready)

    output:
    // Optional: --upload-file writes no report files, --log-file renames the log.
    path('*bulk_report_svn_*.json'), emit: reports, optional: true
    path('*.log'), emit: run_log, optional: true

    script:
    def args = task.ext.args ?: ''
    def repository_root = "${moduleDir}/../../.."
    """
    # where are endrev, startrev, svn defined?
    # assuming we do not need to run this anymore
    # python $params.rfamprod/scripts/apicuron/extract_svn_info.py --end-rev \$endrev --start-rev \$startrev --svn \$svn

    # python $params.rfamprod/scripts/apicuron/bulk_upload.py -f bulk_report_svn.json
    # startrev: svn log -r "{2023-08-17}:HEAD" --limit 1
    # do we need to run this in the datamover queue?
    export PYTHONPATH=${repository_root}
    python3 ${repository_root}/scripts/apicuron/bulk_upload.py ${args}
    """
}

workflow apicuron {
    take:
    start

    main:
    APICURON_REPORT(start)

    emit:
    reports = APICURON_REPORT.out.reports
    run_log = APICURON_REPORT.out.run_log
}

workflow {
    if (params.apicuron) {
        apicuron(channel.of('start'))
    }
}
