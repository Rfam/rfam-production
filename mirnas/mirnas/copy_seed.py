import argparse
import os
import shutil
import subprocess

from scripts.mirnas.mirna_config import NEW_DIR, SEARCH_DIRS, UPDATE_DIR
from scripts.mirnas.update_mirnas_helpers import get_mirna_dict, get_mirna_ids


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
    if new:
        check = create_new_dir(mirna_dir)
        checkout_dir = os.path.join(NEW_DIR, mirna_dir)
    else:
        rfam_acc = mirnas_dict[mirna].keys()[0]
        check = checkout_family(rfam_acc)
        checkout_dir = os.path.join(UPDATE_DIR, rfam_acc)
    if check:
        if mirna.find("relabelled") == -1:
            mirna_dir = mirna + "_relabelled"
        for search_dir in SEARCH_DIRS:
            if os.path.exists(os.path.join(search_dir, mirna_dir)):
                updated_mirna_loc = os.path.join(search_dir, mirna_dir)
                updated_seed = os.path.join(updated_mirna_loc, "SEED")
                if os.path.exists(os.path.join(checkout_dir, "SEED")):
                    os.rename(
                        os.path.join(checkout_dir, "SEED"),
                        os.path.join(checkout_dir, "SEED_old"),
                    )
                shutil.copyfile(updated_seed, os.path.join(checkout_dir, "SEED"))
            else:
                continue


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Checkout (for update) or create dir (for new) for a family, "
        "and copy over the SEED file"
    )
    parser.add_argument(
        "--input",
        help="TSV file with miRNA ID, and threshold value of families to update, "
        "file will also include Rfam acc number if families to update",
    )
    parser.add_argument(
        "--new",
        help="True if miRNA IDs are new families",
        action="store_true",
        default=False,
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()

    if not args.input:
        print("Please provide an input.")
    elif args.new:
        mirna_ids = get_mirna_ids(args.input)
        for mirna_id in mirna_ids:
            copy_seed_file(mirna_id, new=True)
    else:
        mirnas_dict = get_mirna_dict(args.input)
        for mirna_id in mirnas_dict:
            copy_seed_file(mirna_id)
