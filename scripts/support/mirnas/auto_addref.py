import os
import argparse
import subprocess
from scripts.support.mirnas.update_mirnas_helpers import UPDATE_DIR, MEMORY, CPU, LSF_GROUP, get_rfam_accs


def auto_addref(rfam_accs, reference):
    cmd = ("bsub -M {mem} -o {out_file} -e {err_file} -n {cpu} -g {lsf_group} -q production-rh74 "
           "-J {job_name} \"cd {family_dir} && add_ref.pl {ref}\"")
    for family in rfam_accs:
        family_dir = os.path.join(UPDATE_DIR, family)
        if os.path.exists(family_dir):
            lsf_err_file = os.path.join(family_dir, "auto_add_ref.err")
            lsf_out_file = os.path.join(family_dir, "auto_add_ref.out")
            job_name = family
            subprocess.call(
                cmd.format(
                    mem=MEMORY, out_file=lsf_out_file, err_file=lsf_err_file, cpu=CPU, lsf_group=LSF_GROUP,
                    job_name=job_name, family_dir=family_dir, ref=reference), shell=True)
        else:
            continue


def add_ref_sequentially(rfam_accs, reference):
    cmd = ("add_ref.pl {0}".format(reference))

    for family in rfam_accs:
        family_dir = os.path.join(UPDATE_DIR, family)
        if os.path.exists(family_dir):
            os.chdir(family_dir)
            subprocess.call(cmd, shell=True)


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref", help="A string indicating the PubMed id to use for reference",
                        action="store", default="30423142")
    parser.add_argument("--sequential", help="Modify DESC files sequentially",
                        action="store_true", default=False)

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()
    rfam_accs = get_rfam_accs()
    if args.sequential:
        add_ref_sequentially(rfam_accs, args.ref)
    else:
        auto_addref(rfam_accs, args.ref)
