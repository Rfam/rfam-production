import os
import shutil

from scripts.mirnas.mirna_config import COPY_DIR, STK_DIRS


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
    Copy the .stk file to the new miRNA family search directory
    :param mirna: MiRBase ID
    """
    check = create_new_dir(mirna)
    new_search_dir = os.path.join(COPY_DIR, mirna)
    if check:
        stk_file = "{mirna_id}.stk".format(mirna_id=mirna)
        for stk_dir in STK_DIRS:
            stk_file = os.path.join(stk_dir, stk_file)
            if os.path.exists(stk_file):
                # create file named SEED and copy .stk contents
                with open(os.path.join(new_search_dir, "SEED"), "w") as new_seed_file:
                    pass
                shutil.copyfile(stk_file, os.path.join(new_search_dir, "SEED"))
            else:
                continue


if __name__ == "__main__":
    mirna_ids = []
    for directory in STK_DIRS:
        filenames = next(os.walk(directory), (None, None, []))[2]
        for file_name in filenames:
            mirna_id = file_name.replace(".stk", "")
            mirna_ids.append(mirna_id)
    for mirna_id in mirna_ids:
        copy_file(mirna_id)
