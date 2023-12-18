"""
This scripts extracts information from the Rfam SVN repository, and then uploads to APICURON, using the APICURON API. 
Example usage:
    bulk_upload.py --end 16099 --start 15938 --svn /path/to/svn
"""

import argparse
import logging
import requests
import datetime
import json
import re
import subprocess

from config import rfam_config
import scripts.apicuron.conf as conf


def get_author(revision, svn_url):
    """
    Run svnlook to get the author of the commit
    :param revision: svn repo revision number
    :param svn_url: filepath of the svn
    :return: author
    """

    author_cmd = "svnlook author -r {rev} {url}".format(rev=revision, url=svn_url)
    output = subprocess.check_output(author_cmd, shell=True, text=True)
    return output.strip()


def get_family(revision, svn_url):
    """
    Run svnlook to get the directory changed in the commit, then parse this to extract the family name
    :param revision: svn repo revision number
    :param svn_url: filepath of the svn
    :return: Rfam family name
    """

    family_cmd = "svnlook dirs-changed -r {rev} {url}".format(rev=revision, url=svn_url)
    output = subprocess.check_output(family_cmd, shell=True, text=True)
    if match := re.search(r"(RF\d{5})", output):
        return "https://rfam.org/family/" + match.group(1)
    return ""


def get_term(revision, svn_url):
    """
    Run svnlook to get the message associated with the commit
    :param revision: svn repo revision number
    :param svn_url: filepath of the svn
    :return: the activity term, e.g. 'create_family'
    """

    message_cmd = "svnlook log -r {rev} {url}".format(rev=revision, url=svn_url)
    output = subprocess.check_output(message_cmd, shell=True, text=True)
    activity_term = ""
    for ci_term, apicuron_term in conf.checkin_terms.items():
        if ci_term in output:
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
    output = subprocess.check_output(date_cmd, shell=True, text=True)
    timestamp = output[:19]
    return timestamp


def write_report(args):
    """Write the bulk JSON report of activity entries
    Two reports will be written, one which will contain entries for the curators in the
    `svn_authors_current` dict in `conf.py`, and the other will contain entries for the curators in the
    `svn_authors_other` dict in `conf.py`.

    TODO: create only one file

    :param args: CLI arguments of svn url, revision start and revision end
    """

    reports_current = []
    reports_other = []
    url = args.svn

    for rev in range(args.start, args.end):
        author = get_author(rev, url)
        author_orcid = conf.curator_orcids.get(author, author)
        entry = {
            "activity_term": get_term(rev, url),
            "timestamp": get_timestamp(rev, url),
            "curator_orcid": author_orcid,
            "entity_uri": get_family(rev, url),
        }
        if any(value == "" for value in entry.values()):
            continue
        else:
            if author in conf.curator_orcids:
                reports_current.append(entry)
            else:
                reports_other.append(entry)

    today_date = str(datetime.date.today())
    reports_file = "bulk_report_svn_" + today_date + ".json"

    with open(reports_file, "w") as bulk_report:
        reports = {"resource_id": "rfam", "reports": reports_current}
        json.dump(reports, bulk_report, indent=4, sort_keys=True)

    with open("others_" + reports_file, "w") as bulk_report:
        reports = {"resource_id": "rfam", "reports": reports_other}
        json.dump(reports, bulk_report, indent=4, sort_keys=True)

    return reports


def get_header():
    """
    Return header information.
    """

    token = rfam_config.APICURON_TOKEN
    if token is None:
        raise RuntimeError("No APICURON_TOKEN in config")
    header = {
        "version": "2",
        "authorization": token,
    }
    return header


def upload_bulk_report(report):
    """
    Upload a bulk report to APICURON.
    """
    
    try:
        response = requests.post(
            rfam_config.BULK_REPORT_URL,
            files={"reports": open(report, "r")},
            headers=get_header(),
            data={"delete_all": "rfam"},
        )
        response.raise_for_status()
    except requests.exceptions.HTTPError as e:
        logging.debug("HTTP error has occurred uploading to APICURON")
        raise e


def parse_args():
    """
    Parse the CLI arguments
    """
    
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--end", type=int, help="most recent revision number", action="store"
    )
    parser.add_argument(
        "--start",
        type=int,
        help="revision number to start from e.g. revision at last release",
        action="store",
    )
    parser.add_argument(
        "--svn", type=str, help="path to SVN repo to query", action="store"
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    report = write_report(args)
    upload_bulk_report(report)
