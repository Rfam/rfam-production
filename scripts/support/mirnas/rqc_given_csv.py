import os
import argparse
import subprocess

import time

from scripts.support.mirnas.update_mirnas_helpers import (get_rfam_accs, MEMORY, CPU,
                                                          LSF_GROUP)
from scripts.support.mirnas.mirna_config import UPDATE_DIR

families_with_seed_error = []
ignore_seed = []
passed = []
did_not_pass = []


def num_lines_in_seed_file(acc):
    seed_file = os.path.join(UPDATE_DIR, acc, 'SEED')
    with open(seed_file) as f:
        num_lines = sum(1 for _ in f)
    print("Number of lines in SEED file: {0}".format(num_lines))
    return num_lines


def check_seed_match(num_errors, family):
    num_seed_sequences = num_lines_in_seed_file(family)
    percentage_seed_error = (float(num_errors) / float(num_seed_sequences) * 100)
    if percentage_seed_error < 25:
        return True
    else:
        return False


def check_rqc_passes(family):
    family_dir = os.path.join(UPDATE_DIR, family)
    lsf_err_file = os.path.join(family_dir, "auto_rqc.err")

    with open(lsf_err_file) as rqc_output:
        read_output = rqc_output.read()
        if "Family passed with no serious errors" in read_output:
            passed.append(rfam_acc)
            continue_to_check_in = True
        elif "in SEED in not in SCORES list!" in read_output:
            families_with_seed_error.append(rfam_acc)
            number_error_occurrences = read_output.count("in SEED in not in SCORES list!")
            if check_seed_match(number_error_occurrences, rfam_acc):
                print(
                    "At least 75% of sequences match, continuing to check in family {0} with -i seed".format(rfam_acc))
                ignore_seed.append(rfam_acc)
                continue_to_check_in = True
            else:
                print("Error - SEED not in SCORES list. Larger than 25% mismatch. "
                      "Manual intervention required, see error output:{0}".format(lsf_err_file))
                continue_to_check_in = False
        else:
            continue_to_check_in = False
            print("Family {0} failed QC checks. Manual Intervention is required. Please see error output: {1}".format(
                family, lsf_err_file))

    return continue_to_check_in


def run_qc_check(family):
    family_dir = os.path.join(UPDATE_DIR, family)
    lsf_err_file = os.path.join(family_dir, "auto_rqc.err")
    lsf_out_file = os.path.join(family_dir, "auto_rqc.out")
    job_name = os.path.basename(family)
    cmd = ("bsub -M {mem} -o {out_file} -e {err_file} -n {cpu} -g {lsf_group} -J {job_name} "
           "\"cd {update_dir} && rqc-all.pl {rfam_acc}\"")
    subprocess.call(
        cmd.format(
            mem=MEMORY, out_file=lsf_out_file, err_file=lsf_err_file, cpu=CPU, lsf_group=LSF_GROUP,
            job_name=job_name, update_dir=UPDATE_DIR, rfam_acc=family), shell=True)


def write_to_files():
    with open(os.path.join(UPDATE_DIR, 'qc_passed.txt'), 'w') as outfile:
        for family in passed:
            outfile.write(family + '\n')
    with open(os.path.join(UPDATE_DIR, 'qc_passed_ignore_seed.txt'), 'w') as outfile:
        for family in ignore_seed:
            outfile.write(family + '\n')
    with open(os.path.join(UPDATE_DIR, 'qc_not_passed.txt'), 'a') as outfile:
        for family in did_not_pass:
            outfile.write(family + '\n')
    with open(os.path.join(UPDATE_DIR, 'seed_not_in_scores.txt'), 'a') as outfile:
        for family in families_with_seed_error:
            outfile.write(family + '\n')


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv-input",
                        help="CSV file with miRNA id, rfam accession number, threshold value of families to update")

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()
    rfam_accs = get_rfam_accs(csv_file=args.csv_input)
    for rfam_acc in rfam_accs:
        run_qc_check(rfam_acc)
        time.sleep(60)
        if check_rqc_passes(rfam_acc):
            print("Continue to check in {0}".format(rfam_acc))
        else:
            did_not_pass.append(rfam_acc)
            print("Cannot check in {0} - manual intervention required.".format(rfam_acc))

    write_to_files()
