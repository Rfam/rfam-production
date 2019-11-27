import os
import sys
import shutil
import subprocess

def clean_up_sequences(updir):
    """

    :param updir:
    :return:
    """

    seq_dir = os.path.join(updir, "sequences")
    contents = os.listdir(seq_dir)
    for item in contents:
        item_path = os.path.join(seq_dir, item)
        if os.path.isfile(item_path):
            os.remove(item_path)
        elif os.path.isdir(item_path):
            shutil.rmtree(item_path)

if __name__ == '__main__':

    updir = sys.argv[1]

    clean_up_sequences(updir)