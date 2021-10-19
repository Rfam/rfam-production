import json
import os
import subprocess
import argparse

from scripts.support.mirnas.update_mirnas_helpers import MEMORY, CPU, LSF_GROUP, get_mirna_dict
from scripts.support.mirnas.config import UPDATE_DIR, SEARCH_DIRS


def rfmake_serial(family_dir, threshold):
    """
    Call rfmake.pl sequentially
    :param family_dir: location of valid family directory
    :param threshold: manually selected threshold
    """
    serial_cmd = "rfmake.pl -t {threshold} -forcethr"
    os.chdir(family_dir)
    subprocess.call(serial_cmd.format(threshold), shell=True)


def rfmake(family_dir, entry_id, threshold):
    """
    Submit LSF jobs with rfmake.pl
    :param family_dir: location of valid family directory
    :param entry_id: ID of the family to update, miRBase ID or Rfam accession number
    :param threshold: manually selected threshold
    """
    cmd = ("bsub -M {mem} -o {out_file} -e {err_file} -n {cpu} -g {lsf_group} -q production-rh74 "
           "-J {job_name} \"cd {family_dir} && rfmake.pl -t {threshold} -forcethr -a\"")
    lsf_err_file = os.path.join(family_dir, "auto_rfmake.err")
    lsf_out_file = os.path.join(family_dir, "auto_rfmake.out")
    subprocess.call(
        cmd.format(
            mem=MEMORY, out_file=lsf_out_file, err_file=lsf_err_file, cpu=CPU, lsf_group=LSF_GROUP,
            job_name=entry_id, family_dir=family_dir, threshold=threshold), shell=True)


def autorfmake(entryids_thresholds, serial=False):
    """
    Call rfmake.pl with given entry ID and threshold value
    :param entryids_thresholds: key:value pairs of rfam_acc:threshold
    :param serial: run rfmake sequentially if True, else submit to LSF
    """
    could_not_update = []
    family_dir = ""

    for entry in entryids_thresholds:
        entry_id = entry.keys()[0]
        threshold = entry.values()[0]
        if entry_id.startswith('RF'):
            family_dir = os.path.join(UPDATE_DIR, entry_id)
        elif entry_id.startswith('MIPF'):
            for searchdir in SEARCH_DIRS:
                if entry_id.find("relabelled") == -1:
                    family_dir = os.path.join(searchdir, entry_id + "_relabelled")
                else:
                    family_dir = os.path.join(searchdir, entry_id)
        if os.path.exists(family_dir):
            if serial is True:
                rfmake_serial(family_dir, threshold)
            else:
                rfmake(family_dir, entry_id, threshold)
        else:
            could_not_update.append(family_dir)
            continue


def parse_arguments():
    parser = argparse.ArgumentParser(description="run rfmake.pl with a manually selected threshold")
    required_arguments = parser.add_argument_group("required arguments")
    required_arguments.add_argument("--csv-input",
                                    help="CSV file with miRNA id, rfam accession number, "
                                         "threshold value of families to update")
    parser.add_argument("--serial", help="Serial execution of rfmake", action="store_true", default=False)
    parser.add_argument("--thresholds", help="A json file with miRNA : threshold pairs", action="store")

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()
    if args.csv_input:
        mirnas_dict = get_mirna_dict(args.input)
        ids_thresholds = mirnas_dict.values()
    elif args.thresholds:
        json_file = args.thresholds
        with open(json_file, 'r') as fp:
            ids_thresholds = json.load(fp)
    autorfmake(ids_thresholds, args.serial)
