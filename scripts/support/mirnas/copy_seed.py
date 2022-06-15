import argparse
import json
import os
import shutil
import subprocess

from scripts.support.mirnas.update_mirnas_helpers import get_mirna_dict, get_mirna_ids
from scripts.support.mirnas.mirna_config import UPDATE_DIR, SEARCH_DIRS, NEW_DIR


def create_new_dir(mirna):
    """
    Create a new directory for the new miRNA family
    :param mirna: miRBase ID
    :return: True if successful
    """
    family_dir = os.path.join(NEW_DIR, mirna)
    os.mkdir(family_dir)
    if os.path.exists(family_dir):
        return True
    return False


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


def copy_seed_file(mirna, new=False):
    """
    Copy the SEED file from the miRNA search directories to the Rfam family directory
    :param new: if new family
    :param mirna: MiRBase ID
    """

    mirna_dir = mirna
    if mirna.find("relabelled") == -1:
        mirna_dir = mirna + "_relabelled"
    if new:
        check = create_new_dir(mirna_dir)
        checkout_dir = os.path.join(NEW_DIR, mirna_dir)
    else:
        rfam_acc = mirnas_dict[mirna].keys()[0]
        check = checkout_family(rfam_acc)
        checkout_dir = os.path.join(UPDATE_DIR, rfam_acc)
    if check:
        for search_dir in SEARCH_DIRS:
            if os.path.exists(os.path.join(search_dir, mirna_dir)):
                updated_mirna_loc = os.path.join(search_dir, mirna_dir)
                updated_seed = os.path.join(updated_mirna_loc, "SEED")
                if os.path.exists(os.path.join(checkout_dir, "SEED")):
                    os.rename(os.path.join(checkout_dir, "SEED"), os.path.join(checkout_dir, "SEED_old"))
                shutil.copyfile(updated_seed, os.path.join(checkout_dir, "SEED"))
            else:
                continue


def parse_arguments():
    parser = argparse.ArgumentParser(description="Checkout (for update) or create dir (for new) for a family, "
                                                 "and copy over the SEED file")
    mutually_exclusive = parser.add_mutually_exclusive_group()
    mutually_exclusive.add_argument("--csv-input",
                                    help="CSV file with miRNA id, rfam accession number, "
                                         "threshold value of families to update")
    mutually_exclusive.add_argument("--input", help="JSON file with miRNA id : threshold value pairs")

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()
    if args.csv_input:
        mirnas_dict = get_mirna_dict(args.csv_input)
        for mirna_id in mirnas_dict:
            copy_seed_file(mirna_id)
    elif args.input:
        mirna_ids = get_mirna_ids(args.input)
        for mirna_id in mirna_ids:
            copy_seed_file(mirna_id, new=True)
    else:
        print("Please provide an input.")
