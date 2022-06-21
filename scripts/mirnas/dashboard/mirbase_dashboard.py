"""
Generate a dashboard summarising miRBase search results that can be uploaded
to Google Sheets for interactive analysis.

Usage:
python mirbase_dashboard.py <document_id> <sheet_id>
"""

import argparse
import csv
import os

from compute_dashboard_entries import run_rfmake, get_overlaps, is_inconsistent_ss, \
    get_id_matches, get_action
from dashboard_exceptions import SpeciesFileNotFound
from format_dashboard import get_header_line, get_input_data, format_rfam_overlaps, \
    format_mirbase_url, format_rfam_ids, get_summary, get_google_sheets_data
from microrna_progress import updated_families, new_commits
from getters import get_output_path, get_output_url, get_family_location, get_report_url, \
    get_mirbase_alignments, get_mirbase_id, get_mirbase_acc

HTML_REPORTS = '/nfs/public/rw/xfam/rfam/' + 'searches/mirbase'
OUTPUT_FILENAME = 'mirbase-dashboard.tsv'
BLACK_LIST = [
    'MIPF0000338__mir-680'
    'MIPF0000419__mir-574',
    'MIPF0000901__mir-3470',
    'MIPF0001768__mir-7407',  # rfmake running
]  # giant families that cause crashes


def qc_check(new, updated):
    """
    Ensure that all families that are recorded as updated have the correct
    status in the dashboard.

    Found 972 new families but expected to find 979
    Found 267 updated families but expected to find 272

    (Pdb) set(new_commits) - new
    set(['RF03770', 'RF04088', 'RF03311'])
    (Pdb) set(updated_families) - updated
    set(['RF00663', 'RF00859', 'RF00762', 'RF00727', 'RF00716'])
    """
    new_commits_set = set(new_commits)
    updated_families_set = set(updated_families)
    if len(new) != len(new_commits_set):
        msg = 'Found {} new families but expected to find {}'.format(len(new), len(new_commits_set))
        print(new_commits_set - new)
    else:
        msg = 'Found expected number of new families'
    print(msg)
    if len(updated) != len(updated_families_set):
        msg = 'Found {} updated families but expected to find {}'.format(len(updated), len(updated_families_set))
        print(updated_families_set - updated)
    else:
        msg = 'Found expected number of updated families'
    print(msg)


def generate_dashboard(f_out, data, nocache):
    """
    """
    csvwriter = csv.writer(f_out, delimiter='\t', quoting=csv.QUOTE_ALL)
    csvwriter.writerow(get_header_line())
    done_new = set()
    done_updated = set()
    for row_id, identifier in enumerate(get_mirbase_alignments()):
        action = ''
        overlaps = []
        print(identifier)
        if identifier in BLACK_LIST:
            print('Warning: family {} is skipped because it will cause a crash'.format(identifier))
            skip = True
        else:
            skip = False

        if identifier not in data:
            action = 'Fix missing data'
            metadata = {
                'score': '',
                'author': '',
                'comment': '',
            }
        else:
            metadata = data[identifier]
        score = metadata['score']
        location = get_family_location(identifier)
        if is_inconsistent_ss(location):
            action = 'Inconsistent SS_CONS'
            skip = True
        try:
            if not skip:
                run_rfmake(location, score)
        except SpeciesFileNotFound as e:
            print(e.message)
            action = 'Fix rfmake problems'
        try:
            if not action and not skip:
                overlaps = get_overlaps(identifier, nocache)
        except:
            overlaps = []
            action = 'Fix overlap problems'
        overlaps_by_id = get_id_matches(get_mirbase_id(identifier))
        if not overlaps and overlaps_by_id:
            rfam_matches = format_rfam_overlaps(overlaps_by_id)
            rfam_matches = rfam_matches.replace('")', ' Match by ID, no overlap")')
        else:
            rfam_matches = format_rfam_overlaps(overlaps)
        if not action:
            action = get_action(identifier, location, overlaps, overlaps_by_id, score)
        fields = [
            get_summary(row_id, 'label'),
            get_summary(row_id, 'formula'),
            get_report_url(identifier),
            str(score),
            metadata['author'],
            action,
            format_mirbase_url(get_mirbase_acc(identifier)),
            get_mirbase_id(identifier),
            format_rfam_ids(overlaps) if overlaps else format_rfam_ids(overlaps_by_id),
            rfam_matches,
            metadata['comment'],
        ]
        csvwriter.writerow(fields)
        if action == 'Done (new family)':
            done_new.add(overlaps[0] if overlaps else overlaps_by_id[0])
        elif action == 'Done (updated family)':
            done_updated.add(overlaps[0] if overlaps else overlaps_by_id[0])
    qc_check(done_new, done_updated)


def main():
    """
    Main entry point.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('google_doc_id', type=str, help='Google Doc ID', action='store')
    parser.add_argument('google_sheet_id', type=str, help='Google Sheet ID', action='store')
    parser.add_argument("--nocache", help="Recompute cached overlap files", action="store_true", default=False)
    args = parser.parse_args()
    google_doc_id = args.google_doc_id
    google_sheet_id = args.google_sheet_id
    nocache = args.nocache

    google_file = get_google_sheets_data(google_doc_id, google_sheet_id)
    data = get_input_data(google_file)
    with open(get_output_path(), 'w') as f_out:
        generate_dashboard(f_out, data, nocache)
    os.system('rm {}'.format(google_file))
    print("""
    Created a file on disk: {path}
    The file is available at: {url}

    To load the data in Google Sheets:

    1. Save the file locally

    wget -O {filename} {url}

    2. Manually import file {filename} into Google Sheets.

    3. Follow instructions in Readme.

    """.format(path=get_output_path(), url=get_output_url(OUTPUT_FILENAME),
               filename=OUTPUT_FILENAME))


if __name__ == '__main__':
    main()
