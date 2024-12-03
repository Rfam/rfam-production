import argparse
import json
import os
import subprocess

from scripts.mirnas.mirna_config import CPU, LSF_GROUP, MEMORY, NEW_DIR, UPDATE_DIR
from scripts.mirnas.update_mirnas_helpers import get_id_thresholds, get_mirna_dict


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
    cmd = (
        "bsub -M {mem} -q short -o {out_file} -e {err_file} -n {cpu} -g {lsf_group} "
        '-J {job_name} "cd {family_dir} && rfmake.pl -t {threshold} -forcethr -a -local"'
    )
    lsf_err_file = os.path.join(family_dir, "auto_rfmake.err")
    lsf_out_file = os.path.join(family_dir, "auto_rfmake.out")
    subprocess.call(
        cmd.format(
            mem=MEMORY,
            out_file=lsf_out_file,
            err_file=lsf_err_file,
            cpu=CPU,
            lsf_group=LSF_GROUP,
            job_name=entry_id,
            family_dir=family_dir,
            threshold=threshold,
        ),
        shell=True,
    )


def autorfmake(entryids_thresholds, serial=False):
    """
    Call rfmake.pl with given entry ID and threshold value
    :param entryids_thresholds: key:value pairs of rfam_acc:threshold
    :param serial: run rfmake sequentially if True, else submit to LSF
    """
    could_not_update = []
    family_dir = ""

    for entry_id, threshold in entryids_thresholds.items():
        if entry_id.startswith("RF"):
            family_dir = os.path.join(UPDATE_DIR, entry_id)
        elif entry_id.startswith("MIPF"):
            family_dir = os.path.join(NEW_DIR, entry_id)
        if os.path.exists(family_dir):
            if serial is True:
                rfmake_serial(family_dir, threshold)
            else:
                rfmake(family_dir, entry_id, threshold)
        else:
            could_not_update.append(family_dir)
            continue
    if could_not_update:
        print("Could not update: {0}".format(could_not_update))


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="run rfmake.pl with a manually selected threshold"
    )
    required_arguments = parser.add_argument_group("required arguments")
    required_arguments.add_argument(
        "--input",
        help="TSV file with miRNA ID, and threshold value of families to update, "
        "file will also include Rfam acc number if families to update",
    )
    parser.add_argument(
        "--thresholds", help="A json file with miRNA : threshold pairs", action="store"
    )
    parser.add_argument(
        "--serial",
        help="Serial execution of rfmake",
        action="store_true",
        default=False,
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
    # new families
    if args.input and args.new:
        ids_thresholds = get_id_thresholds(args.input)
    # update families
    elif args.input:
        mirnas_dict = get_mirna_dict(args.input)
        ids_thresholds = mirnas_dict.values()
        ids_thresholds = {k: v for d in ids_thresholds for k, v in d.items()}
    elif args.thresholds:
        json_file = args.thresholds
        with open(json_file, "r") as fp:
            ids_thresholds = json.load(fp)
    else:
        print("Please provide an input.")

    autorfmake(ids_thresholds, args.serial)
