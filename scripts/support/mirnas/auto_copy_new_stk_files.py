import argparse
import os
import shutil

from scripts.support.mirnas.mirna_config import STK_DIRS

COPY_DIR = "/nfs/production/agb/rfam/RELEASES/14.9/microrna/batch3_chunk2_searches"


def create_new_dir(mirna):
    """
    Create a new directory for the new miRNA family
    :param mirna: miRBase ID
    :return: True if successful
    """
    family_dir = os.path.join(COPY_DIR, mirna)
    os.mkdir(family_dir)
    if os.path.exists(family_dir):
        return True
    return False


def copy_file(mirna):
    """

    :param mirna: MiRBase ID
    """
    check = create_new_dir(mirna)
    checkout_dir = os.path.join(COPY_DIR, mirna)
    if check:
        for new_stk_dir in STK_DIRS:
            if os.path.exists(os.path.join(new_stk_dir, mirna)):
                stk_file = "{mirna_id}.stk".format(mirna_id=mirna)
                updated_mirna_loc = os.path.join(new_stk_dir, mirna)
                updated_seed = os.path.join(updated_mirna_loc, stk_file)
                if os.path.exists(os.path.join(checkout_dir, stk_file)):
                    # rename .stk file SEED
                    os.rename(os.path.join(checkout_dir, stk_file), os.path.join(checkout_dir, "SEED"))
                shutil.copyfile(updated_seed, os.path.join(checkout_dir, "SEED"))
            else:
                continue


def parse_arguments():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--input", help="List of miRNA IDs")

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()
    mirna_ids = []
    for dir in STK_DIRS:
        for file_name in dir:
            id = file_name.replace('.stk', '')
            mirna_ids.append(id)
    for mirna_id in mirna_ids:
        copy_file(mirna_id)
    else:
        print("Please provide an input.")
