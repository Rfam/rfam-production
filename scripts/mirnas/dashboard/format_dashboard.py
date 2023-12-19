import csv
import os
import tempfile
import time
import collections as coll

from getters import get_rfam_id, get_rfam_clan


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
               'Rfam Clan(s)', 'Comment (editable)', ]
    return columns


def get_input_data(filename):
    """
    Read current spreadsheet to get manually curated fields.
    """
    data = {}
    csvreader = csv.reader(open(filename, 'r'), delimiter='\t')
    next(csvreader)  # skip header
    for row in csvreader:
        data[row[2]] = {
            'score': row[3],
            'author': row[4],
            'comment': row[10],
        }

    return data


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


def format_rfam_ids(overlaps):
    """
    Given a list of Rfam accessions, return a list of Rfam IDs.
    """
    rfam_ids = []
    for rfam_acc in sorted(overlaps):
        rfam_ids.append(get_rfam_id(rfam_acc))
    return ', '.join(rfam_ids)


def format_rfam_clans(overlaps):
    """
    Given a list of Rfam accessions, return a string of CLAN_ID(RFAM_IDS).

    >>> rfam_rfam_clans(['RF00818', 'RF00716'])
    'CL00084 (RF00818, RF00716)'
    >>> rfam_rfam_clans(['RF00716', 'RF00070'])
    'CL00051 (RF00070), CL00084 (RF00716)'
    >>> rfam_rfam_clans(['RF00716', 'RF00070', 'RF02913'])
    'CL00051 (RF00070), CL00084 (RF00716), No clan (RF02913)'
    """
    clans = coll.defaultdict(set)
    for rfam_acc in overlaps:
        clan_id = get_rfam_clan(rfam_acc) or 'No clan'
        clans[clan_id].add(str(rfam_acc))

    sorted_clans = []
    for clan, rfam_accs in clans.items():
        rfam_accs = filter(None, rfam_accs)
        sorted_clans.append("{clan}({accs})".format(
                            clan=clan,
                            accs=','.join(sorted(rfam_accs))))
    return ', '.join(sorted(sorted_clans))


def get_summary(row_id, data_type):
    """
    The Google Sheets formulas are used to preserve interactivity of the dashboard.
    It is possible to compute these numbers from the data itself but the dashboard
    would have to be recomputed once the fields are updated.
    """
    summary = [
        {
            'label': 'Total',
            'formula': '=COUNTA(C:C)-1',  # subtract the header row
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


def get_google_sheets_data(document_id, sheet_id):
    """
    Download current dashboard data from Google Sheets and store in a temporary
    file to avoid dealing with permissions.
    """
    temp_file = tempfile.NamedTemporaryFile(delete=False)
    cmd = 'wget -O "{}" "https://docs.google.com/spreadsheets/d/{}/export?format=tsv&gid={}"'
    os.system(cmd.format(temp_file.name, document_id, sheet_id))
    return temp_file.name
