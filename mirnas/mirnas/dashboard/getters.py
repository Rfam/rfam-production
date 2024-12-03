import glob
import os
import random
import string

import requests
from dashboard_config import HTML_REPORTS, OUTPUT_FILENAME, SEARCH_DIR


def get_output_path():
    """
    Get a full path to the output file on disk.
    The path is mapped to a public URL on the preview website.
    """
    return os.path.join(HTML_REPORTS, OUTPUT_FILENAME)


def get_output_url(filename):
    """
    Get a public URL for the output file.
    Include a random cache-busting string to avoid accidental re-download
    of a previously generated file (especially useful when debugging).
    """
    return "/".join(
        [
            "https://preview.rfam.org",
            HTML_REPORTS.split(os.sep)[-2],
            HTML_REPORTS.split(os.sep)[-1],
            filename
            + "?cachebust-"
            + "".join(
                random.choice(string.ascii_uppercase + string.digits) for _ in range(10)
            ),
        ]
    )


def get_family_location(identifier):
    """
    Get a path to the search files for a given identifier, including checking
    for several possible filenames.
    """
    location = os.path.join(SEARCH_DIR, identifier)
    if os.path.exists(location):
        return location
    location = location + "_relabelled"
    if os.path.exists(location):
        return location
    raise Exception("Search location not found")


def get_report_url(identifier):
    """
    Get a URL of the report on the preview website.
    """
    file_path = get_report_path(identifier)
    if file_path:
        url = file_path.replace(
            HTML_REPORTS, "https://preview.rfam.org/searches/mirbase"
        )
        return '=HYPERLINK("{}","{}")'.format(url, identifier)
    else:
        return identifier


def get_report_path(identifier):
    """
    Check if an html report exists for a given search identifier.
    """
    report = os.path.join(HTML_REPORTS, identifier + ".html")
    if os.path.exists(report):
        return report
    report = report.replace(".html", "_relabelled.html")
    if os.path.exists(report):
        return report
    return None


def get_rfam_id(rfam_acc):
    """
    Get Rfam ID for an Rfam accession.
    """
    url = "http://www.ebi.ac.uk/ebisearch/ws/rest/rfam?query={}%20AND%20entry_type:%22Family%22&fields=name&format=json"
    data = requests.get(url.format(rfam_acc))
    try:
        if data.json()["hitCount"] != 0:
            return data.json()["entries"][0]["fields"]["name"][0]
        else:
            return ""
    except KeyError as e:
        print(
            "Warning: response from {url} did not contain a hit count, {e}".format(
                url=url.format(rfam_acc), e=e
            )
        )
        return ""


def get_rfam_clan(rfam_acc):
    """
    Get Rfam ID for an Rfam accession.
    """
    url = "http://www.ebi.ac.uk/ebisearch/ws/rest/rfam?query={}%20AND%20entry_type:%22Clan%22&fields=name&format=json"
    data = requests.get(url.format(rfam_acc))
    if data.json()["hitCount"] != 0:
        return data.json()["entries"][0]["fields"]["name"][0]
    return None


def get_mirbase_alignments():
    """
    Get a list of alignments provided by miRBase.
    """
    script_location = os.path.dirname(os.path.abspath(__file__))
    data_location = os.path.join(script_location, "mirbase-seeds", "*.stk")
    files = glob.glob(data_location)
    return [
        os.path.basename(x.replace(".stk", "").replace("_relabelled", ""))
        for x in files
    ]


def get_mirbase_id(identifier):
    """
    Example:
    MIPF0000024__mir-103

    >>> get_mirbase_id('MIPF0000024__mir-103')
    'mir-103'
    >>> get_mirbase_id('MIPF0000024')
    ''
    >>> get_mirbase_id('MIPF0000024_mir-103')
    ''
    >>> get_mirbase_id('MIPF0000024__mir_103')
    'mir_103'
    """
    if "__" in identifier:
        _, mirbase_id = identifier.split("__")
    else:
        mirbase_id = ""
    return mirbase_id


def get_mirbase_acc(identifier):
    """
    Example:
    MIPF0000024__mir-103

    >>> get_mirbase_id('MIPF0000024__mir-103')
    'MIPF0000024'
    >>> get_mirbase_id('MIPF0000024')
    ''
    >>> get_mirbase_id('MIPF0000024_mir-103')
    ''
    >>> get_mirbase_id('MIPF0000024__mir_103')
    'MIPF0000024'
    """
    if "__" in identifier:
        mirbase_acc, _ = identifier.split("__")
    else:
        mirbase_acc = ""
    return mirbase_acc
