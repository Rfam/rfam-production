import datetime
import json
import re
import subprocess

import argparse

import scripts.apicuron.conf as conf


def get_author(revision, svn_url):
    """
    Run svnlook to get the author of the commit
    :param revision: svn repo revision number
    :param svn_url: filepath of the svn
    :return: author
    """

    author_cmd = "svnlook author -r {rev} {url}".format(rev=revision, url=svn_url)
    p = subprocess.Popen(author_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output, stderr = p.communicate()
    return output.strip()


def get_family(revision, svn_url):
    """
    Run svnlook to get the directory changed in the commit, then parse this to extract the family name
    :param revision: svn repo revision number
    :param svn_url: filepath of the svn
    :return: Rfam family name
    """

    family_cmd = "svnlook dirs-changed -r {rev} {url}".format(rev=revision, url=svn_url)
    p = subprocess.Popen(family_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output, stderr = p.communicate()
    match = re.search(r'RF\d{5}', output)
    family = ''
    if match:
        family = match.group(0)
    return family


def get_term(revision, svn_url):
    """
    Run svnlook to get the message associated with the commit
    :param revision: svn repo revision number
    :param svn_url: filepath of the svn
    :return: the activity term, e.g. 'create_family'
    """

    message_cmd = "svnlook log -r {rev} {url}".format(rev=revision, url=svn_url)
    p = subprocess.Popen(message_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output, stderr = p.communicate()
    activity_term = ''
    for rfci_term, apicuron_term in conf.checkin_terms.items():
        if rfci_term in output:
            activity_term = apicuron_term
    return activity_term


def get_timestamp(revision, svn_url):
    """
    Run svnlook to get the timestamp of the commit
    :param revision: svn repo revision number
    :param svn_url: filepath of the svn
    :return: timestamp, as string
    """

    date_cmd = "svnlook date -r {rev} {url}".format(rev=revision, url=svn_url)
    p = subprocess.Popen(date_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output, stderr = p.communicate()
    timestamp = output[:19]
    return timestamp


def parse_args():
    """
    Parse the CLI arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--end', type=int, help='most recent revision number', action='store')
    parser.add_argument('--start', type=int, help='revision number to start from e.g. revision at last release',
                        action='store')
    parser.add_argument('--svn', type=str, help='SVN repo to query', action='store')
    return parser.parse_args()


def main():
    """
    Example usage:
    extract_svn_info.py --end 16099 --start 15938 --svn /path/to/svn
    """
    args = parse_args()
    reports_current = []
    reports_other = []
    url = args.svn
    for rev in range(args.start_rev, args.end_rev):
        author = get_author(rev, url)
        if any(author == a for a in conf.svn_authors_current):
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
                'curator_orcid': conf.curator_orcids[author] if conf.curator_orcids[author] else author,
                'entity_uri': "https://rfam.org/family/" + get_family(rev, url)
            }
            reports_other.append(entry)
    today_date = str(datetime.date.today())
    reports_file = 'bulk_report_svn_' + today_date + '.json'
    with open(reports_file, 'w') as bulk_report:
        reports = {'resource_id': 'rfam', 'reports': reports_current}
        json.dump(reports, bulk_report, indent=4, sort_keys=True)
    with open('others_' + reports_file, 'w') as bulk_report:
        reports = {'resource_id': 'rfam', 'reports': reports_other}
        json.dump(reports, bulk_report, indent=4, sort_keys=True)


if __name__ == '__main__':
    main()
