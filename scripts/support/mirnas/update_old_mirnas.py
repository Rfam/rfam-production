import argparse
import os
import shutil
import subprocess

from scripts.support.mirnas.update_mirnas_helpers import UPDATE_DIR, get_data_from_csv

searchdirs = ["/hps/nobackup/production/xfam/rfam/RELEASES/14.3/miRNA_relabelled/batch1_chunk1_searches",
              "/hps/nobackup/production/xfam/rfam/RELEASES/14.3/miRNA_relabelled/batch1_chunk2_searches",
              "/hps/nobackup/production/xfam/rfam/RELEASES/14.3/miRNA_relabelled/batch2/searches"]


def checkout_family(rfam_acc):
    os.chdir(UPDATE_DIR)
    cmd = "rfco.pl {0}".format(rfam_acc)
    subprocess.call(cmd, shell=True)
    if os.path.exists(os.path.join(UPDATE_DIR, rfam_acc)):
        return True
    return False


def copy_seed_file(mirna):
    mirna_dir = mirna
    if mirna.find("relabelled") == -1:
        mirna_dir = mirna + "_relabelled"
    rfam_acc = mirnas_dict[mirna].keys()
    check = checkout_family(rfam_acc)
    checkout_dir = os.path.join(UPDATE_DIR, rfam_acc)
    if check:
        for search_dir in searchdirs:
            if os.path.exists(os.path.join(search_dir, mirna_dir)):
                updated_mirna_loc = os.path.join(search_dir, mirna_dir)
                updated_seed = os.path.join(updated_mirna_loc, "SEED")
                os.rename(os.path.join(checkout_dir, "SEED"), os.path.join(checkout_dir, "SEED_old"))
                shutil.copyfile(updated_seed, os.path.join(checkout_dir, "SEED"))
            else:
                continue


def parse_arguments():
    parser = argparse.ArgumentParser(description='Script to update miRNAs')
    required_arguments = parser.add_argument_group("required arguments")
    # input is now the csv file
    required_arguments.add_argument("--input", help="CSV file", type=str)

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()
    mirnas_dict = get_data_from_csv(args.input)
    for mirna_id in mirnas_dict:
        copy_seed_file(mirna_id)
