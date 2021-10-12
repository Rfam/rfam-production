import argparse
import fileinput
import os
import sys

from scripts.support.mirnas.update_mirnas_helpers import get_rfam_accs, UPDATE_DIR

field_options = {
    'AU': 'AU   Griffiths-Jones SR; 0000-0001-6043-807X\n',
    'SE': 'SE   Griffiths-Jones SR\n',
    'SS': 'SS   Predicted; RNAalifold\n'
}


def replace_field(field):
    for line in fileinput.input('DESC', inplace=True):
        if line.strip().startswith(field):
            line = field_options[field]
        sys.stdout.write(line)


def update_desc_fields(rfam_accessions):
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
    required_arguments.add_argument("--input", help="CSV file", type=str)

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()
    rfam_accs = get_rfam_accs(args.input)
    update_desc_fields(rfam_accs)
