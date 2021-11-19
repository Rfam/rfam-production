import argparse
from datetime import datetime
import fileinput

HEADER = \
    "!version date: {date} \n" \
    "!description: A mapping of GO terms to Rfam release <XX.X>.\n" \
    "!external resource: https://rfam.org/\n" \
    "!citation: Kalvari et al. (2020) Nucl. Acids Res. 49: D192-D200\n" \
    "!contact: rfam-help@ebi.ac.uk\n" \
    "!\n"


def add_header(rfam2go_file):
    """
    Add header to rfam2go with version date
    """
    header = HEADER.format(date=datetime.today().strftime('%Y-%m-%d'))

    fi = fileinput.input(rfam2go_file, inplace=1)
    for line in fi:
        if fi.isfirstline():
            print header + line,
        else:
            print line,


def parse_arguments():
    parser = argparse.ArgumentParser()
    required_arguments = parser.add_argument_group("required arguments")
    required_arguments.add_argument("--input", help="")

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()
    add_header(args.input)
