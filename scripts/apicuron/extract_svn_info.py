import json
import re
import subprocess

import argparse

import scripts.apicuron.conf as conf


def get_author(revision, svn_url):
    author_cmd = "svnlook author -r {rev} {url}".format(rev=revision, url=svn_url)
    p = subprocess.Popen(author_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output, stderr = p.communicate()
    return output.strip()


def get_family(revision, svn_url):
    family_cmd = "svnlook dirs-changed -r {rev} {url}".format(rev=revision, url=svn_url)
    p = subprocess.Popen(family_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output, stderr = p.communicate()
    match = re.search(r'RF\d{5}', output)
    family = ''
    if match:
        family = match.group(0)
    return family


def get_term(revision, svn_url):
    message_cmd = "svnlook log -r {rev} {url}".format(rev=revision, url=svn_url)
    p = subprocess.Popen(message_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output, stderr = p.communicate()
    activity_term = ''
    for rfci_term, apicuron_term in conf.rfci_terms.items():
        if rfci_term in output:
            activity_term = apicuron_term
    return activity_term


def get_timestamp(revision, svn_url):
    date_cmd = "svnlook date -r {rev} {url}".format(rev=revision, url=svn_url)
    p = subprocess.Popen(date_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output, stderr = p.communicate()
    timestamp = output[:19]
    return timestamp


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--end-rev', type=int, help='most recent revision number', action='store')
    parser.add_argument('--start-rev', type=int, help='revision number to start from e.g. revision at last release',
                        action='store')
    parser.add_argument('--svn', type=str, help='SVN repo to query', action='store')
    args = parser.parse_args()
    reports_current = []
    reports_other = []
    url = args.svn
    for rev in range(args.start_rev, args.end_rev):
        author = get_author(rev, url)
        if any(author == a for a in conf.svn_authors):
            entry = {
                'activity_term': get_term(rev, url),
                'timestamp': get_timestamp(rev, url),
                'curator_orcid': conf.curator_orcids[author],
                'entity_uri': "https://rfam.org/family/" + get_family(rev, url)
            }
            reports_current.append(entry)
        else:
            entry = {
                'activity_term': get_term(rev, url),
                'timestamp': get_timestamp(rev, url),
                'curator_orcid': author,
                'entity_uri': "https://rfam.org/family/" + get_family(rev, url)
            }
            reports_other.append(entry)
    with open('bulk_report_svn.json', 'w') as bulk_report:
        reports = {'resource_id': 'rfam', 'reports': reports_current}
        json.dump(reports, bulk_report, indent=4, sort_keys=True)
    with open('bulk_report_svn_others.json', 'w') as bulk_report:
        reports = {'resource_id': 'rfam', 'reports': reports_other}
        json.dump(reports, bulk_report, indent=4, sort_keys=True)


if __name__ == '__main__':
    main()
