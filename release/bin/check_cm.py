import argparse
import sys
from pathlib import Path


def num_families_check(num_accs: int, statfile: Path) -> bool:
    """
    Check that Rfam.cm contains the correct number of families
    :param num_accs: number of families
    :param statfile: output of cmstat
    """
    count = 0

    with statfile.open("r") as stats:
        for line in stats:
            if "#" not in line:
                count += 1

    print("Number of families: {0}".format(count))
    if num_accs != count:
        print("Number of families in CM file is incorrect!")
        return False
    return True


def desc_check(num_accs: int, cmfile: Path) -> bool:
    """
    Check the number of DESC lines - should be 2 * number of families
    :param num_accs: number of families
    :param cmfile: CM file
    """
    count = 0

    with cmfile.open("r") as cm:
        for line in cm:
            if "DESC" in line:
                count += 1

    print("Number of DESC lines: {0}".format(count))
    if (2 * num_accs) != count:
        print("Number of DESC lines is incorrect!")
        return False
    return True


def all_accs_in_cm(accs: ty.List[str], cmfile: Path) -> bool:
    """
    Check that all accessions in the database are present in the Rfam.cm file before we continue
    :param accs: list of all rfam_acc found in the family table of the rfam_live database
    :param cm: Rfam.cm file
    :return:
    """

    missing = set()
    with cmfile.open("r") as cm:
        contents = cm.read()
        for acc in accs:
            if acc not in contents:
                missing.add(acc)
    if missing:
        missed = ", ".join(missing)
        print("CM file missing: {}".format(missed))
        return False
    return True


def parse_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-ac", "--accessions", help="Accession file", required=True)
    parser.add_argument("-cf", "--cm-file", help="CM file", required=True)
    parser.add_argument("-sf", "--stat-file", help="cmstat output", required=True)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    accessions = args.accessions
    success = num_families_check(len(accessions), args.stat_file)
    success &= desc_check(len(accessions), args.cm_file)
    success &= all_accs_in_cm(accessions, args.cm_file)
    if not success:
        print("At least one check failed")
        sys.exit(1)
