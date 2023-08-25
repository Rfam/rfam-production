from dataclasses import dataclass
import typing as ty
import datetime
import json
import re
import subprocess

import argparse

import scripts.apicuron.conf as conf


@dataclass
class RfamActivity:
    acitivity_term: str
    timestamp: str
    curator_orcid: ty.Optional[str]
    entity_uri: str


def get_author(revision: int, svn_url: str) -> str:
    """
    Run svnlook to get the author of the commit
    :param revision: svn repo revision number
    :param svn_url: filepath of the svn
    :return: author
    """

    author_cmd = "svnlook author -r {rev} {url}".format(rev=revision, url=svn_url)
    output = subprocess.check_output(author_cmd, shell=True, text=True)
    return output.strip()


def get_family(revision: int, svn_url: str):
    """
    Run svnlook to get the directory changed in the commit, then parse this to extract the family name
    :param revision: svn repo revision number
    :param svn_url: filepath of the svn
    :return: Rfam family name
    """

    family_cmd = "svnlook dirs-changed -r {rev} {url}".format(rev=revision, url=svn_url)
    output = subprocess.check_output(family_cmd, shell=True, text=True)
    if match := re.search(r"(RF\d{5})", output):
        return match.group(1)
    return ""


def get_term(revision: int, svn_url: str) -> str:
    """
    Run svnlook to get the message associated with the commit
    :param revision: svn repo revision number
    :param svn_url: filepath of the svn
    :return: the activity term, e.g. 'create_family'
    """

    message_cmd = "svnlook log -r {rev} {url}".format(rev=revision, url=svn_url)
    output = subprocess.check_output(message_cmd, shell=True, text=True)
    activity_term = ""
    for rfci_term, apicuron_term in conf.checkin_terms.items():
        if rfci_term in output:
            activity_term = apicuron_term
    return activity_term


def get_timestamp(revision: int, svn_url: str) -> str:
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


def main():
    """
    Example usage:
    extract_svn_info.py --end 16099 --start 15938 --svn /path/to/svn
    """
    args = parse_args()
    reports_current = []
    reports_other = []
    url = args.svn
    for rev in range(args.start, args.end):
        author = get_author(rev, url)
        author_orcid = conf.curator_orcids.get(author, author)
        entry = RfamActivity(
            get_term(rev, url),
            get_timestamp(rev, url),
            author_orcid,
            "https://rfam.org/family/" + get_family(rev, url)
        )

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


if __name__ == "__main__":
    main()
