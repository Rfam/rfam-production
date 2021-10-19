import argparse
import os
import shutil
import subprocess

from scripts.support.mirnas.update_mirnas_helpers import get_mirna_dict
from scripts.support.mirnas.config import UPDATE_DIR, SEARCH_DIRS


def checkout_family(rfam_acc):
    """
    Run rfco.pl to checkout the given Rfam family
    :param rfam_acc: accession number of the family to check out
    :return: True if the new family directory exists after rfco.pl is ran, else False
    """

    os.chdir(UPDATE_DIR)
    cmd = "rfco.pl {0}".format(rfam_acc)
    subprocess.call(cmd, shell=True)
    if os.path.exists(os.path.join(UPDATE_DIR, rfam_acc)):
        return True
    return False


def copy_seed_file(mirna):
    """
    Copy the SEEd file from the miRNA search directories to the Rfam family directory
    :param mirna: MiRBase ID
    """

    mirna_dir = mirna
    if mirna.find("relabelled") == -1:
        mirna_dir = mirna + "_relabelled"
    rfam_acc = mirnas_dict[mirna].keys()[0]
    check = checkout_family(rfam_acc)
    checkout_dir = os.path.join(UPDATE_DIR, rfam_acc)
    if check:
        for search_dir in SEARCH_DIRS:
            if os.path.exists(os.path.join(search_dir, mirna_dir)):
                updated_mirna_loc = os.path.join(search_dir, mirna_dir)
                updated_seed = os.path.join(updated_mirna_loc, "SEED")
                os.rename(os.path.join(checkout_dir, "SEED"), os.path.join(checkout_dir, "SEED_old"))
                shutil.copyfile(updated_seed, os.path.join(checkout_dir, "SEED"))
            else:
                continue


def parse_arguments():
    parser = argparse.ArgumentParser(description="Checkout a family, and copy over the SEED file")
    required_arguments = parser.add_argument_group("required arguments")
    required_arguments.add_argument("--csv-input",
                                    help="CSV file with miRNA id, rfam accession number, "
                                         "threshold value of families to update")

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()
    mirnas_dict = get_mirna_dict(args.csv_input)
    for mirna_id in mirnas_dict:
        copy_seed_file(mirna_id)
