import os
import sys
import argparse
import subprocess
import json

import time
from subprocess import Popen, PIPE

from scripts.support.mirnas.update_mirnas_helpers import get_data_from_csv, get_rfam_accs, UPDATE_DIR, MEMORY, CPU, \
    LSF_GROUP

search_dirs = ["/hps/nobackup/production/xfam/rfam/RELEASES/14.3/miRNA_relabelled/batch1_chunk1_searches",
               "/hps/nobackup/production/xfam/rfam/RELEASES/14.3/miRNA_relabelled/batch1_chunk2_searches",
               "/hps/nobackup/production/xfam/rfam/RELEASES/14.3/miRNA_relabelled/batch2/searches"]


def check_desc_ga(DESC, cut_ga):
    process = Popen(['grep', "GA", DESC], stdin=PIPE, stdout=PIPE, stderr=PIPE)
    output, err = process.communicate()

    if output.find("%.2f" % float(cut_ga)) == -1:
        return False

    return True


def check_family_passes_qc(family_dir):
    dir_elements = os.path.split(family_dir)
    search_dir = dir_elements[0]

    os.chdir(search_dir)

    process = Popen(["rqc-all.pl", dir_elements[1]], stdin=PIPE, stdout=PIPE, stderr=PIPE)
    output = process.communicate()[1]

    if output.find("Family passed with no serious errors") == -1:
        return False

    return True


def find_rqc_error(rqc_output):
    error_types = {"CM": 1, "FORMAT": 1, "OVERLAP": 1, "STRUCTURE": 1,
                   "MISSING": 1, "SEQUENCE": 1, "NON-CODING": 1}

    rqc_lines = [x.strip() for x in rqc_output.split("\n") if x != '']

    error_type_count = len(error_types.keys())

    success_string = "--%s check completed with no major errors"

    for error_type in error_types.keys():
        for rqc_line in rqc_lines:
            if rqc_line.find(success_string % error_type) != -1:
                error_types[error_type] = 0
                break

    return error_types


def fetch_rqc_output(family_dir):
    dir_elements = os.path.split(family_dir)
    parent_dir = dir_elements[0]

    os.chdir(parent_dir)

    process = Popen(["rqc-all.pl", dir_elements[1]], stdin=PIPE, stdout=PIPE, stderr=PIPE)
    output = process.communicate()[1]

    return output


def generate_rqc_report(rqc_error_types):
    # develop functionality here
    pass


def fetch_format_error_lines(rqc_output):
    rqc_lines = [x.strip() for x in rqc_output.split("\n") if x != '']

    error_lines = []

    start = '(2) FORMAT CHECK'
    end = '(3) OVERLAP CHECK'
    flag = False

    for line in rqc_lines:
        if line == start:
            flag = True
        if line == end:
            flag = False
        if flag is True:
            error_lines.append(line)

    return error_lines


def check_error_is_fatal(error_lines):
    error_string = "\n".join(error_lines)

    if error_string.find("FATAL") != -1:
        return True

    return False


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


def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument("--input", help="A CSV file with miRNAs to commit",
                        action="store")
    parser.add_argument("--log-dir", help="Log file destination",
                        action="store", default=os.getcwd())
    parser.add_argument("--overlaps", help="Displays family overlaps",
                        action="store_true", default=False)
    parser.add_argument("--format", help="Displays families failing FORMAT check",
                        action="store_true", default=False)

    return parser.parse_args()


def check_rqc_passes(rfam_acc):
    family_dir = os.path.join(UPDATE_DIR, rfam_acc)
    lsf_out_file = os.path.join(family_dir, "auto_rqc.out")

    with open(lsf_out_file) as rqc_output:
        if "Family passed with no serious errors" in rqc_output:
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


if __name__ == '__main__':

    args = parse_arguments()
    rfam_accs = get_rfam_accs(csv_file=args.input)
    for rfam_acc in rfam_accs:
        run_qc_check(rfam_acc)
        time.sleep(10)
        if check_rqc_passes(rfam_acc):
            print('{0} passed QC checks'.format(rfam_acc))
        else:
            print('{0} DID NOT PASS QC checks. Please check the output at {1}'.format(
                rfam_acc, os.path.join(UPDATE_DIR, rfam_acc, "auto_rqc.out")))

    #             rqc_output = fetch_rqc_output(family_dir_loc)
    #             qc_error_dict = find_rqc_error(rqc_output)
    #             # count = 0
    #             # for key in qc_error_dict.keys():
    #             #	count = count + qc_error_dict[key]
    #             # if count != 0:
    #             # if qc_error_dict["STRUCTURE"] == 1 or qc_error_dict["FORMAT"] == 1:
    #             #	print (family_dir_loc)
    #             if args.overlaps:
    #                 if qc_error_dict["OVERLAP"] == 1:
    #                     overlaps = extract_family_overlaps(rqc_output)
    #                     for family in overlaps:
    #                         print("%s\t%s" % (mirna_d, family))
    #             elif args.format:
    #                 if qc_error_dict["FORMAT"] == 1:
    #                     print(family_dir_loc)

    # print find_2_seed_family(rqc_output)
