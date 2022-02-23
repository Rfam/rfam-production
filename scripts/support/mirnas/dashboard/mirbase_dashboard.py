"""
Generate a dashboard summarising miRBase search results that can be uploaded
to Google Sheets for interactive analysis.

Usage:
python mirbase_dashboard.py <document_id> <sheet_id>
"""

import argparse
import csv
import glob
import os
import re
import tempfile
import time

import requests


from microrna_progress import updated_families, new_commits


SEARCH_DIR = '/hps/nobackup/production/xfam/rfam/RELEASES/14.8/microrna/searches'
HTML_REPORTS = '/nfs/public/rw/xfam/rfam/test/' + 'searches/mirbase'
OUTPUT_FILENAME = 'mirbase-dashboard.tsv'
BLACK_LIST = [
    'MIPF0000419__mir-574',
    'MIPF0000901__mir-3470',
    'MIPF0001768__mir-7407', #rfmake running
] # giant families that cause crashes


def get_output_path():
    """
    Get a full path to the output file on disk.
    The path is mapped to a public URL on the preview website.
    """
    return os.path.join(HTML_REPORTS, OUTPUT_FILENAME)


def get_output_url():
    """
    Get a public URL for the output file.
    """
    return '/'.join([
        'https://preview.rfam.org',
        HTML_REPORTS.split(os.sep)[-2],
        HTML_REPORTS.split(os.sep)[-1],
        OUTPUT_FILENAME
    ])


def get_family_location(identifier):
    """
    Get a path to the search files for a given identifier, including checking
    for several possible filenames.
    """
    location = os.path.join(SEARCH_DIR, identifier)
    if os.path.exists(location):
        return location
    location = location + '_relabelled'
    if os.path.exists(location):
        return location
    raise Exception('Search location not found')


class SpeciesFileNotFound(Exception):
    """
    Raised when rfmake fails.
    """
    pass


def run_rfmake(location, score):
    """
    Run rfmake with a manually selected threshold if not done already.
    """
    if not score or '?' in score:
        return
    score = float(score)
    species_file = os.path.join(location, 'species')
    if not os.path.exists(species_file):
        raise SpeciesFileNotFound('Species file not found in {}'.format(location))
    ga_threshold = None
    with open(species_file, 'r') as f_in:
        for line in f_in:
            match = re.search(r'CURRENT GA THRESHOLD: (.+) BITS', line)
            if match:
                ga_threshold = float(match.group(1))
                break
    if ga_threshold and ga_threshold == score:
        pass
    else:
        cmd = 'cd {} && rfmake.pl -t {} && cd - > /dev/null'.format(location, score)
        os.system(cmd)


def get_overlaps(identifier, nocache):
    """
    Run rqc-overlaps and get overlapping families.
    External overlap [SS] of CM001599.2/64499845-64499913 with RF00600:CM001599.2/64499850-64499995 by 64
    """
    location = get_family_location(identifier)
    overlap_file = os.path.join(location, 'overlap')
    rfam_accs = set()
        cmd = 'cd {} && touch SEED CM DESC TBLOUT SCORES && cd - > /dev/null'.format(location)
        os.system(cmd)
    if not os.path.exists(overlap_file) or os.stat(overlap_file).st_size == 0 or nocache:
        cmd = 'cd {} && rqc-overlap.pl {} && cd - > /dev/null'.format(SEARCH_DIR, os.path.basename(location))
        os.system(cmd)
    with open(overlap_file, 'r') as f_overlap:
        for line in f_overlap:
            match = re.search(r'with (RF\d{5})', line)
            if match:
                rfam_accs.add(match.group(1))
    return list(rfam_accs)


def get_report_url(identifier):
    """
    Get a URL of the report on the preview website.
    """
    file_path = get_report_path(identifier)
    if file_path:
        url = file_path.replace(HTML_REPORTS, 'https://preview.rfam.org/searches/mirbase')
        return '=HYPERLINK("{}","{}")'.format(url, identifier)
    else:
        return identifier


def get_report_path(identifier):
    """
    Check if an html report exists for a given search identifier.
    """
    report = os.path.join(HTML_REPORTS, identifier + '.html')
    if os.path.exists(report):
        return report
    report = report.replace('.html', '_relabelled.html')
    if os.path.exists(report):
        return report
    return None


def is_single_seed(location):
    """
    Check if the seed contains only 1 sequence.
    """
    seed_file = os.path.join(location, 'SEED')
    alistat_file = seed_file + '_alistat.txt'
    cmd = 'esl-alistat {} > {}'.format(seed_file, alistat_file)
    os.system(cmd)
    with open(alistat_file, 'r') as f_alistat:
        for line in f_alistat:
            if line.startswith('Number of sequences:'):
                fields = re.split(r'\:\s+', line)
                if len(fields) == 2:
                    num_seed = int(fields[1])
                    if num_seed == 1:
                        return True
                    else:
                        return False
    return None


def is_inconsistent_ss(location):
    """
    Check if seed alignment has an inconsistent secondary structure.
    """
    seed_file = os.path.join(location, 'SEED')
    alistat_file = seed_file + '_alistat.txt'
    cmd = 'esl-alistat --bpinfo test.bp.txt {} > {} 2>&1'.format(seed_file, alistat_file)
    os.system(cmd)
    with open(alistat_file, 'r') as f_alistat:
        for line in f_alistat:
            if line.startswith('Consensus structure string is inconsistent'):
                return True
    return False


def get_new_or_updated(overlaps):
    """
    Check if this is one of the microRNA family that has already been updated or
    newly committed.
    """
    if len(overlaps) == 1:
        rfam_acc = overlaps[0]
        if rfam_acc in updated_families:
            return 'Updated'
        elif rfam_acc in new_commits:
            return 'New'
    return 'False'


def get_action(identifier, location, overlaps, overlaps_by_id, score):
    """
    Identify what should be done with a family.
    """
    action = ''
    if overlaps:
        overlap_status = get_new_or_updated(overlaps)
    elif overlaps_by_id:
        overlap_status = get_new_or_updated(overlaps_by_id)
    else:
        overlap_status = None

    mirbase_id = get_mirbase_id(identifier)
    if len(overlaps_by_id) == 1 and len(overlaps) == 1 and overlaps[0] != overlaps_by_id[0]:
        return 'Fix ID mismatch'
    elif len(overlaps_by_id) == 1 and len(overlaps) == 1 and overlaps[0] == overlaps_by_id[0]:
        rfam_id = get_rfam_id(overlaps[0])
    elif len(overlaps) == 1 and not overlaps_by_id:
        rfam_id = get_rfam_id(overlaps[0])
    elif len(overlaps_by_id) == 1 and not overlaps:
        rfam_id = get_rfam_id(overlaps_by_id[0])

    if 'HYPERLINK' not in get_report_url(identifier):
        action = 'Generate report'
    elif is_single_seed(location):
        action = '1_SEED'
    elif not score or score == '?':
        action = 'Choose threshold'
    elif overlap_status == 'New' and rfam_id.lower() == mirbase_id.lower():
        action = 'Done (new family)'
    elif overlap_status == 'Updated' and rfam_id.lower() == mirbase_id.lower():
        action = 'Done (updated family)'
    elif score > 0 and not overlaps and not overlaps_by_id:
        action = 'New family'
    elif score > 0 and (overlaps or overlaps_by_id):
        action = 'Update seed'
    return action


def format_rfam_overlaps(overlaps):
    """
    Format Rfam accessions as website links (except when there are several
    accessions as it's not possible to have multiple links in one cell).
    """
    if not overlaps:
        return ''
    if len(overlaps) == 1:
        rfam_acc = overlaps[0]
        return '=HYPERLINK("https://rfam.org/family/{0}", "{0}")'.format(rfam_acc)
    else:
        return ', '.join(list(overlaps))


def format_mirbase_url(identifier):
    """
    Format miRBase hyperlinks.
    """
    link = '=HYPERLINK("http://www.mirbase.org/summary.shtml?fam={0}","{0}")'
    return link.format(identifier)


def get_rfam_id(rfam_acc):
    """
    Get Rfam ID for an Rfam accession.
    """
    url = 'http://www.ebi.ac.uk/ebisearch/ws/rest/rfam?query={}%20AND%20entry_type:%22Family%22&fields=name&format=json'
    data = requests.get(url.format(rfam_acc))
    if data.json()['hitCount'] != 0:
        return(data.json()['entries'][0]['fields']['name'][0])
    else:
        return ''


def format_rfam_ids(overlaps):
    """
    Given a list of Rfam accessions, return a list of Rfam IDs.
    """
    rfam_ids = []
    for rfam_acc in sorted(overlaps):
        rfam_ids.append(get_rfam_id(rfam_acc))
    return ', '.join(rfam_ids)


def get_id_matches(mirbase_id):
    """
    In case there are no overlaps, check for ID matches.
    """
    url = 'http://www.ebi.ac.uk/ebisearch/ws/rest/rfam?query="{}"%20AND%20entry_type:%22Family%22&fields=name&format=json'
    data = requests.get(url.format(mirbase_id))
    if data.json()['hitCount'] == 1:
        rfam_id = data.json()['entries'][0]['id']
        return [rfam_id]
    return []


def format_output_line(fields):
    """
    Format tabular output.
    Skip two columns to give space for the summary formulas.
    """
    return '\t\t' + '\t'.join(fields) + '\n'


def get_header_line():
    """
    Return spreadsheet column headers.
    """
    columns = ['Summary', '', 'Report', 'Score (editable)', 'Author (editable)',
               'Action', 'miRBase AC', 'miRBase ID', 'Rfam ID(s)', 'Rfam AC(s)',
               'Comment (editable)',]
    return columns


def get_input_data(filename):
    """
    Read current spreadsheet to get manually curated fields.
    """
    data = {}
    csvreader = csv.reader(open(filename, 'r'), delimiter='\t')
    next(csvreader) # skip header
    for row in csvreader:
        data[row[2]] = {
            'score': row[3],
            'author': row[4],
            'comment': row[10],
        }

    return data


def get_summary(row_id, data_type):
    """
    The Google Sheets formulas are used to preserve interactivity of the dashboard.
    It is possible to compute these numbers from the data itself but the dashboard
    would have to be recomputed once the fields are updated.
    """
    summary = [
        {
            'label': 'Total',
            'formula': '=COUNTA(C:C)-1', # subtract the header row
        },
        {
            'label': 'Update seed',
            'formula': '=COUNTIF(F:F, "Update seed")',
        },
        {
            'label': 'New family',
            'formula': '=COUNTIF(F:F, "New family")',
        },
        {
            'label': 'Generate report',
            'formula': '=COUNTIF(F:F, "Generate report")',
        },
        {
            'label': '1_SEED',
            'formula': '=COUNTIF(F:F, "1_SEED")',
        },
        {
            'label': 'Choose threshold',
            'formula': '=COUNTIF(F:F, "Choose threshold")',
        },
        {
            'label': 'Inconsistent SS_CONS',
            'formula': '=COUNTIF(F:F, "Inconsistent SS_CONS")',
        },
        {
            'label': 'Fix a problem',
            'formula': '=COUNTIF(F:F, "Fix*")',
        },
        {
            'label': 'Done',
            'formula': '=COUNTIF(F:F, "Done*")',
        },
        {
            'label': 'Anton',
            'formula': '=COUNTIF(E:E, "Anton")',
        },
        {
            'label': 'Ioanna',
            'formula': '=COUNTIF(E:E, "Ioanna")',
        },
        {
            'label': 'Nancy',
            'formula': '=COUNTIF(E:E, "Nancy")',
        },
        {
            'label': 'Sam',
            'formula': '=COUNTIF(E:E, "Sam")',
        },
        {
            'label': 'Blake',
            'formula': '=COUNTIF(E:E, "Blake")',
        },
        {
            'label': 'Emma',
            'formula': '=COUNTIF(E:E, "Emma")',
        },
        {
            'label': 'Last updated',
            'formula': '{}'.format(time.strftime("%d/%m/%Y %H:%M")),
        },
    ]
    if data_type not in ['label', 'formula']:
        return ''
    if row_id < len(summary):
        return summary[row_id][data_type]
    return ''


def get_mirbase_alignments():
    """
    Get a list of alignments provided by miRBase.
    """
    script_location = os.path.dirname(os.path.abspath(__file__))
    data_location = os.path.join(script_location, 'mirbase-seeds', '*.stk')
    files = glob.glob(data_location)
    return [os.path.basename(x.replace('.stk', '').replace('_relabelled', '')) for x in files]


def get_google_sheets_data(document_id, sheet_id):
    """
    Download current dashboard data from Google Sheets and store in a temporary
    file to avoid dealing with permissions.
    """
    temp_file = tempfile.NamedTemporaryFile(delete=False)
    cmd = 'wget -O "{}" "https://docs.google.com/spreadsheets/d/{}/export?format=tsv&gid={}"'
    os.system(cmd.format(temp_file.name, document_id, sheet_id))
    return temp_file.name


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


def get_mirbase_id(identifier):
    """
    Example:
    MIPF0000024__mir-103
    """
    if '__' in identifier:
        _, mirbase_id = identifier.split('__')
    else:
        mirbase_id = ''
    return mirbase_id


def get_mirbase_acc(identifier):
    """
    Example:
    MIPF0000024__mir-103
    """
    if '__' in identifier:
        mirbase_acc, _ = identifier.split('__')
    else:
        mirbase_acc = ''
    return mirbase_acc


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

    """.format(path=get_output_path(), url=get_output_url(),
               filename=OUTPUT_FILENAME))


if __name__ == '__main__':
    main()
