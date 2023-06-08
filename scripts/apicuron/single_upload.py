import logging

import requests

from config import rfam_local as config

# examples
SINGLE_REPORT_PAYLOAD = {
    "activity_term": "update_family",
    "resource_id": "rfam",
    "timestamp": "2023-04-24T09:59:17.000Z",
    "curator_orcid": "0000-0003-2958-6837",
    "entity_uri": "https://identifiers.org/rfam:RF00657"
}

SINGLE_REPORT_URL = 'https://dev.apicuron.org/api/reports/single'


def get_header():
    """
    Return header information.
    """
    token = config.APICURON_TOKEN
    if token is None:
        raise RuntimeError('No APICURON_TOKEN in config')
    header = {
        'accept: application/json',
        'Content-Type: application/json',
        'Authorization: ' + token
    }
    return header


def upload_single_report():
    """
    Upload a single report entry to APICURON.
    """
    try:
        response = requests.post(SINGLE_REPORT_URL, json=SINGLE_REPORT_PAYLOAD, headers=get_header())
        response.raise_for_status()
    except requests.exceptions.HTTPError as e:
        logging.debug('HTTP error has occurred uploading to APICURON')
        raise e
