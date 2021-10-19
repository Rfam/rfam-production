import argparse
import fileinput
import os
import sys

from scripts.support.mirnas.update_mirnas_helpers import get_rfam_accs
from scripts.support.mirnas.config import UPDATE_DIR

field_options = {
    'AU': 'AU   Griffiths-Jones SR; 0000-0001-6043-807X\n',
    'SE': 'SE   Griffiths-Jones SR\n',
    'SS': 'SS   Predicted; RNAalifold\n'
}


def replace_field(field):
    """
    Replace the line in the DESC file for the given field
    :param field: AU, SE, or SS - field line to update
    """
    for line in fileinput.input('DESC', inplace=True):
        if line.strip().startswith(field):
            line = field_options[field]
        sys.stdout.write(line)


def update_desc_fields(rfam_accessions):
    """
    Update the DESC fields AU, SEE, and SS
    :param rfam_accessions:
    :return:
    """
    fields = ['AU', 'SE', 'SS']
    for family in rfam_accessions:
        family_dir = os.path.join(UPDATE_DIR, family)
        if os.path.exists(family_dir):
            os.chdir(family_dir)
            for field in fields:
                replace_field(field)


def parse_arguments():
    parser = argparse.ArgumentParser()
    required_arguments = parser.add_argument_group("required arguments")
    required_arguments.add_argument(
        "--csv-input", help="CSV file with miRNA id, rfam accession number, threshold value of families to update")

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()
    rfam_accs = get_rfam_accs(args.input)
    update_desc_fields(rfam_accs)
