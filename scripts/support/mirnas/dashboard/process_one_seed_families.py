
"""
Generate an archive with files for families with a single seed (1_SEED).

Usage:
python process_one_seed_family.py <document_id> <sheet_id>
"""


import argparse
import csv
import os
import re
import tempfile


from mirbase_dashboard import get_family_location, get_google_sheets_data, \
                              get_output_url, HTML_REPORTS
from find_family_overlaps import parse_outlist_file


OUTPUT_FILENAME = 'one_seed.tar.gz'


def get_input_data(filename):
    """
    Parse input file and keep only those families that are marked
    as `Update seed`.
    """
    accessions = {}
    with open(filename, 'r') as fp:
        csvreader = csv.reader(fp, delimiter='\t', quoting=csv.QUOTE_ALL)
        for fields in csvreader:
            if fields[5] != '1_SEED':
                continue
            id_text = fields[2]
            if 'HYPERLINK' in id_text:
                match = re.search(r',"(.+)"\)', id_text)
                if match:
                    identifier = match.group(1)
                else:
                    identifier = id_text
            else:
                identifier = id_text
            rfam_accs = re.findall(r'RF\d{5}', fields[9])
            accessions[identifier] = {
                'threshold': [fields[3]],
                'rfam_acc': list(set(rfam_accs)),
            }
    return accessions


def process_one_seed_family(output_folder, accessions, identifier):
    """
    """
    threshold = accessions[identifier]['threshold'][0]
    data_path = get_family_location(identifier)
    align = os.path.join(data_path, 'align-{}'.format(threshold))
    align_with_seed = os.path.join(data_path, 'align-with-seed-{}'.format(threshold))
    align_with_seed_pfam = os.path.join(data_path, 'align-with-seed-pfam-{}'.format(threshold))
    species = os.path.join(data_path, 'species')
    outlist = os.path.join(data_path, 'outlist')
    outlist_info = parse_outlist_file(outlist)

    full_hits_found = False
    for _, data in outlist_info.iteritems():
        if data['seq_label'] != 'SEED':
            full_hits_found = True
            break
    if not full_hits_found:
        print('Only seed hits found in {}'.format(identifier))
        return

    if not os.path.exists(align):
        cmd = 'cd {} && rfmake.pl -t {} -local -a -forcethr -relax && cp align {} && cd - > /dev/null'
        os.system(cmd.format(data_path, threshold, os.path.basename(align)))
    if not os.path.exists(align_with_seed) or os.stat(align_with_seed).st_size == 0:
        cmd = 'cd {} && esl-reformat fasta {} | cmalign --mapali SEED CM - > {} && cd - > /dev/null'
        os.system(cmd.format(data_path, os.path.basename(align), os.path.basename(align_with_seed)))
    if not os.path.exists(align_with_seed_pfam) or os.stat(align_with_seed_pfam).st_size == 0:
        cmd = 'esl-reformat pfam {} > {}'.format(align_with_seed, align_with_seed_pfam)
        os.system(cmd)
    output_family_folder = os.path.join(output_folder, identifier)
    cmd = 'mkdir -p {}'.format(output_family_folder)
    os.system(cmd)
    for filename in [species, outlist, align_with_seed_pfam]:
        cmd = 'cp {} {}'.format(filename, output_family_folder)
        os.system(cmd)


def main():
    """
    Main entry point.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('google_doc_id', type=str, help='Google Doc ID', action='store')
    parser.add_argument('google_sheet_id', type=str, help='Google Sheet ID', action='store')
    args = parser.parse_args()
    google_doc_id = args.google_doc_id
    google_sheet_id = args.google_sheet_id

    google_file = get_google_sheets_data(google_doc_id, google_sheet_id)
    accessions = get_input_data(google_file)

    output_folder = tempfile.mkdtemp()
    print(output_folder)
    for _, identifier in enumerate(accessions.keys()):
        if 'MIPF0001768__mir-7407' in identifier:
            print('Only seed hits found in {}'.format(identifier))
            continue
        process_one_seed_family(output_folder, accessions, identifier)

    filename = os.path.join(HTML_REPORTS, OUTPUT_FILENAME)
    cmd = 'tar -cvzf {} {}'.format(filename, output_folder)
    os.system(cmd)
    os.system('rm {}'.format(google_file))
    os.system('rm -Rf {}'.format(output_folder))
    print("""
    Created a file on disk: {path}
    The file is available at: {url}
    """.format(path=filename, url=get_output_url(OUTPUT_FILENAME)))


if __name__ == "__main__":
    main()
