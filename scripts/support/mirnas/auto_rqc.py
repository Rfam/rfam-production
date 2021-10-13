import os
import argparse
import subprocess

import time

from scripts.support.mirnas.update_mirnas_helpers import (get_rfam_accs, UPDATE_DIR, MEMORY, CPU,
                                                          LSF_GROUP)


def extract_family_overlaps(rqc_output):
    error_lines = []

    rqc_lines = [x.strip() for x in rqc_output.split("\n") if x != '']

    start = '(3) OVERLAP CHECK'
    end = '(4) STRUCTURE CHECK'
    flag = False

    # read error lines
    for line in rqc_lines:
        if line == start:
            flag = True
        if line == end:
            flag = False
        if flag is True:
            error_lines.append(line)

    overlap_accessions = {}
    # process error lines
    for line in error_lines:
        rfam_acc = ""
        if line.find("RF") != -1:
            rfam_acc = line.split(':')[0][-7:]
            if rfam_acc not in overlap_accessions:
                overlap_accessions[rfam_acc] = ''

    return overlap_accessions.keys()


def check_rqc_passes(rfam_acc):
    family_dir = os.path.join(UPDATE_DIR, rfam_acc)
    lsf_err_file = os.path.join(family_dir, "auto_rqc.err")

    with open(lsf_err_file) as rqc_output:
        read_output = rqc_output.read()
        if "Family passed with no serious errors" in read_output:
            return True

    return False


def run_qc_check(rfam_acc):
    family_dir = os.path.join(UPDATE_DIR, rfam_acc)
    lsf_err_file = os.path.join(family_dir, "auto_rqc.err")
    lsf_out_file = os.path.join(family_dir, "auto_rqc.out")
    job_name = os.path.basename(rfam_acc)
    cmd = ("bsub -M {mem} -o {out_file} -e {err_file} -n {cpu} -g {lsf_group} -J {job_name} "
           "\"cd {update_dir} && rqc-all.pl {rfam_acc}\"")
    subprocess.call(
        cmd.format(
            mem=MEMORY, out_file=lsf_out_file, err_file=lsf_err_file, cpu=CPU, lsf_group=LSF_GROUP,
            job_name=job_name, update_dir=UPDATE_DIR, rfam_acc=rfam_acc), shell=True)


def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument("--input", help="A CSV file*",
                        action="store")
    parser.add_argument("--overlaps", help="Displays family overlaps",
                        action="store_true", default=False)

    return parser.parse_args()


if __name__ == '__main__':

    args = parse_arguments()
    rfam_accs = get_rfam_accs(csv_file=args.input)
    did_not_pass = []
    passed = []
    for rfam_acc in rfam_accs:
        run_qc_check(rfam_acc)
        time.sleep(120)
        if check_rqc_passes(rfam_acc):
            print('{0} passed QC checks'.format(rfam_acc))
            passed.append(rfam_acc)
        else:
            did_not_pass.append(rfam_acc)
            print('{0} DID NOT PASS QC checks. Please check the output at {1}'.format(
                rfam_acc, os.path.join(UPDATE_DIR, rfam_acc, "auto_rqc.err")))
    print("Families that did not pass QC checks: {0}".format(did_not_pass))
    with open(os.path.join(UPDATE_DIR, 'qc_passed.txt'), 'w') as outfile:
        for family in passed:
            outfile.write(family)
    with open(os.path.join(UPDATE_DIR, 'did_not_pass_qc.txt'), 'w') as outfile:
        for family in did_not_pass:
            outfile.write(family)
