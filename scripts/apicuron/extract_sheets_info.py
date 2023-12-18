"""Extract information from Google Sheets

This scripts extracts information from the miRNA summary dashboard, which is a Google Sheets file. The script was initially 
required to prepare the info from this sheet for APICURON submissions. It is not required to run regularly, as we now
have a regaulr way to upload APICURON info via bulk_upload.py. However, if any information is lost or needs to be 
recreated, this script will be useful. 
"""

import argparse
import json

import pandas as pd
from faker import Faker

from scripts.apicuron.conf import curator_orcids, sheet_terms

doc_url = "https://docs.google.com/spreadsheets/d/{doc_id}/gviz/tq?tqx=out:csv&sheet={sheet_id}"


def get_random_timestamp():
    """
    Get a random timestamp since the start of the year
    :return: random timestamp 
    """
    fake = Faker()
    random_timestamp = fake.date_time_between(start_date="-2y", end_date="now")
    return random_timestamp.strftime("%Y-%m-%dT%H:%M:%SZ")


def get_report_entries(doc_id, sheet_id):
    """
    Parse the Google sheets doc to get the entries to add to the reports file
    :param doc_id: ID of the Google sheet
    :param sheet_id: ID of the individual sheet
    :return: dict of the entries
    """
    url = doc_url.format(doc_id=doc_id, sheet_id=sheet_id)
    df = pd.read_csv(url)
    reports = []
    for index, row in df.iterrows():
        entry = {
            "activity_term": sheet_terms[row["Action"]],
            "timestamp": get_random_timestamp(),
            "curator_orcid": curator_orcids[row["Author"]],
            "entity_uri": "https://rfam.org/family/" + row["Rfam AC(s)"][:7],
        }
        reports.append(entry)
    return reports


def parse_args():
    """
    Parse the CLI arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("doc_id", type=str, help="Google Doc ID", action="store")
    parser.add_argument("sheet_id", type=str, help="Google Sheet ID", action="store")
    return parser.parse_args()


def main():
    args = parse_args()
    reports = get_report_entries(args.doc_id, args.sheet_id)
    with open("bulk_report_from_sheets.json", "w") as report:
        reports = {"resource_id": "rfam", "reports": reports}
        json.dump(reports, report, indent=4, sort_keys=True)


if __name__ == "__main__":
    main()
