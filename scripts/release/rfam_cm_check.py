import argparse
import sys

from scripts.export.generate_ftp_files import get_all_rfam_accessions


def num_families_check(num_accs, statfile):
    """
    Check that Rfam.cm contains the correct number of families
    :param num_accs: number of families
    :param statfile: output of cmstat
    """
    count = 0

    with open(statfile) as f:
        for line in f:
            if '#' not in line:
                count += 1
    print("Number of families: {0}".format(count))

    if num_accs != count:
        print("Number of families in CM file is incorrect!")
        sys.exit()


def desc_check(num_accs, cmfile):
    """
    Check the number of DESC lines - should be 2 * number of families
    :param num_accs: number of families
    :param cmfile: CM file
    """
    count = 0

    with open(cmfile) as f:
        for line in f:
            if 'DESC' in line:
                count += 1
    print("Number of DESC lines: {0}".format(count))

    if (2*num_accs) != count:
        print("Number of DESC lines is incorrect!")
        sys.exit()


def parse_args():
    """
    Parse the cli arguments when calling this script to insert a text file to the PDB table in the database.
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-cf', '--cm-file', help='CM file', required=True)
    parser.add_argument('-sf', '--stat-file', help='cmstat output', required=True)
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    accessions = get_all_rfam_accessions()
    num_families_check(len(accessions), args.stat_file)
    desc_check(len(accessions), args.cm_file)
