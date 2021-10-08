import os
import subprocess
import argparse

from scripts.support.mirnas.update_mirnas_helpers import UPDATE_DIR, MEMORY, CPU, LSF_GROUP, get_data_from_csv


def autorfmake(accessions_thresholds, serial=False):
    cmd = ("bsub -M {mem} -o {out_file} -e {err_file} -n {cpu} -g {lsf_group} -q production-rh74 "
           "-J {job_name} \"cd {family_dir} && rfmake.pl -t {threshold} -forcethr -a\"")
    serial_cmd = "rfmake.pl -t {threshold} -forcethr"
    could_not_update = []
    for entry in accessions_thresholds:
        rfam_acc = entry.keys()[0]
        threshold = entry.keys()[1]
        family_dir = os.path.join(UPDATE_DIR, rfam_acc)
        if os.path.exists(family_dir):
            if serial is True:
                os.chdir(family_dir)
                subprocess.call(serial_cmd.format(threshold), shell=True)
            else:
                lsf_err_file = os.path.join(family_dir, "auto_rfmake.err")
                lsf_out_file = os.path.join(family_dir, "auto_rfmake.out")
                subprocess.call(
                    cmd.format(
                        mem=MEMORY, out_file=lsf_out_file, err_file=lsf_err_file, cpu=CPU, lsf_group=LSF_GROUP,
                        job_name=rfam_acc, family_dir=family_dir, threshold=threshold), shell=True)
        else:
            could_not_update.append(family_dir)
            continue


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--serial", help="Serial execution of rfmake", action="store_true", default=False)

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()
    mirnas_dict = get_data_from_csv()
    # key:value pairs of rfam_acc:threshold
    ids_thresholds = mirnas_dict.values()
    autorfmake(ids_thresholds, args.serial)
