import json

import argparse
from faker import Faker
import pandas as pd

from scripts.apicuron.conf import sheet_terms, curator_orcids

doc_url = 'https://docs.google.com/spreadsheets/d/{doc_id}/gviz/tq?tqx=out:csv&sheet={sheet_id}'


def get_random_timestamp():
    fake = Faker()
    random_timestamp = fake.date_time_between(start_date='-1y', end_date='now')
    return random_timestamp.strftime('%Y-%m-%dT%H:%M:%SZ')


def write_report(reports):
    with open('bulk_report.json', 'w') as bulk_report:
        json.dump(reports, bulk_report, indent=4, sort_keys=True)


def get_report_entries(doc_id, sheet_id):

    url = doc_url.format(doc_id=doc_id, sheet_id=sheet_id)
    df = pd.read_csv(url)
    reports = []
    for index, row in df.iterrows():
        entry = {
            'activity_term': sheet_terms[row['Action']],
            'timestamp': get_random_timestamp(),
            'curator_orcid': curator_orcids[row['Author']],
            'entity_uri': 'https://rfam.org/family/' + row['Rfam AC(s)'][:7]
        }
        reports.append(entry)
    return reports


def main():
    """
    Main entry point.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('doc_id', type=str, help='Google Doc ID', action='store')
    parser.add_argument('sheet_id', type=str, help='Google Sheet ID', action='store')
    args = parser.parse_args()
    doc_id = args.google_doc_id
    sheet_id = args.google_sheet_id
    reports = get_report_entries(doc_id, sheet_id)
    write_report(reports)


if __name__ == '__main__':
    main()
