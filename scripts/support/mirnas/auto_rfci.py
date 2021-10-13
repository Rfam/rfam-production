import os
import argparse
import subprocess

from scripts.support.mirnas.update_mirnas_helpers import UPDATE_DIR, MEMORY, CPU, LSF_GROUP


def check_in(rfam_acc):
    family_dir = os.path.join(UPDATE_DIR, rfam_acc)
    lsf_err_file = os.path.join(family_dir, "auto_rfci.err")
    lsf_out_file = os.path.join(family_dir, "auto_rfci.out")
    job_name = rfam_acc
    cmd = ("bsub -M {mem} -o {out_file} -e {err_file} -n {cpu} -g {lsf_group} -J {job_name} "
           "\"cd {update_dir} && rfci.pl -m 'Update using miRBase seed' {rfam_acc}\"")
    subprocess.call(
        cmd.format(
            mem=MEMORY, out_file=lsf_out_file, err_file=lsf_err_file, cpu=CPU, lsf_group=LSF_GROUP,
            job_name=job_name, update_dir=UPDATE_DIR, rfam_acc=rfam_acc), shell=True)


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="list of families to check in", action="store")

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()
    families_to_check_in = args.input
    for family in families_to_check_in:
        check_in(family)

