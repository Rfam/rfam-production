import argparse
import sys

from utils.db_utils import fetch_all_rfam_accs


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

    if (2 * num_accs) != count:
        print("Number of DESC lines is incorrect!")
        sys.exit()


def all_accs_in_cm(accs, cm):
    """
    Check that all accessions in the database are present in the Rfam.cm file before we continue
    :param accs: list of all rfam_acc found in the family table of the rfam_live database
    :param cm: Rfam.cm file
    :return:
    """

    with open(cm, 'r') as f:
        contents = f.read()
        for acc in accs:
            if acc not in contents:
                print("Accession not found in cm file, please ensure this cm file has been created: {0}".format(acc))


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
    accessions = fetch_all_rfam_accs()
    num_families_check(len(accessions), args.stat_file)
    desc_check(len(accessions), args.cm_file)
    all_accs_in_cm(accessions, args.cm_file)
