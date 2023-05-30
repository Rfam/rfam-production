import argparse
import logging

import requests

from config import rfam_config


def get_header():
    """
    Return header information.
    """
    token = rfam_config.APICURON_TOKEN
    if token is None:
        raise RuntimeError('No APICURON_TOKEN in config')
    header = {
        'version': '2',
        'authorization': token,
    }
    return header


def upload_bulk_report(report):
    """
    Upload a bulk report to APICURON.
    """
    try:
        response = requests.post(rfam_config.BULK_REPORT_URL, files={'reports': open(report, 'r')}, headers=get_header())
        response.raise_for_status()
    except requests.exceptions.HTTPError as e:
        logging.debug('HTTP error has occurred uploading to APICURON')
        raise e


def parse_args():
    """
    Parse the CLI arguments
    """
    parser = argparse.ArgumentParser(description='Bulk report submission to apicuron')
    parser.add_argument('-f', '--file', help='JSON file with data to upload', required=True)
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    if args.file:
        upload_bulk_report(args.file)
    else:
        raise Exception('Please provide a JSON file to upload data to APICURON')
