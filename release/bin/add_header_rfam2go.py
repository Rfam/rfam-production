import argparse
import fileinput
from datetime import datetime

HEADER = """
!version date: {date}
!description: A mapping of GO terms to Rfam release {version}
!external resource: https://rfam.org/
!citation: Kalvari et al. (2020) Nucl. Acids Res. 49: D192-D200
!contact: rfam-help@ebi.ac.uk
!
"""


def add_header(rfam2go_file, release_version):
    """
    Add header to rfam2go with version date
    """
    header = HEADER.format(
        date=datetime.today().strftime("%Y-%m-%d"), version=release_version
    )

    print(header)
    fileinput
    fi = fileinput.input(rfam2go_file, inplace=1)
    for line in fi:
        if fi.isfirstline():
            print(header + line)
        else:
            print(line)


def parse_arguments():
    parser = argparse.ArgumentParser()
    required_arguments = parser.add_argument_group("required arguments")
    required_arguments.add_argument("--input", help="rfam2fo file")
    required_arguments.add_argument("--rfam-version", help="release version number")

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()
    add_header(args.input, args.version)
