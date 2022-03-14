import os
import subprocess
import argparse

from scripts.support.mirnas.update_mirnas_helpers import get_rfam_accs
from scripts.support.mirnas.mirna_config import UPDATE_DIR, MEMORY, LSF_GROUP


def launch_new_rfsearch(family_dir, cpu):
    """
    Launches a new LSF job with rfsearch.pl
    :param family_dir: The location of a valid family directory (directory with Rfam ID)
    :param cpu: Number of CPUs to use per thread
    """

    lsf_err_file = os.path.join(family_dir, "auto_rfsearch.err")
    lsf_out_file = os.path.join(family_dir, "auto_rfsearch.out")
    job_name = os.path.basename(family_dir)

    if os.path.exists(os.path.join(family_dir, "DESC")) is False:
        option = '-nodesc'
    else:
        option = '-ignoresm'

    file_extensions = ["cmsearch", "tbl", "err"]

    for file_type in file_extensions:
        subprocess.call("rm -f %s/*.%s" % (family_dir, file_type), shell=True)

    cmd = ("bsub -M {mem} -o {out_file} -e {err_file} -n {cpu} -g {lsf_group} -q production-rh74 "
           "-J {job_name} \"cd {family_dir} && rfsearch.pl -t 25 -cnompi -q production-rh74 -relax {option}\"")
    subprocess.call(
        cmd.format(
            mem=MEMORY, out_file=lsf_out_file, err_file=lsf_err_file, cpu=cpu, lsf_group=LSF_GROUP,
            job_name=job_name, family_dir=family_dir, option=option), shell=True)


def remove_all_gaps(family_dir):
    """
    Rename SEED to OLD_SEED, then regenerate SEED file removing all gaps
    :param family_dir: The location of a valid family directory (directory with Rfam ID)
    :return: True if the gaps are removed and the SEED is updated, else False
    """

    seed = os.path.join(family_dir, "SEED")
    old_seed = os.path.join(family_dir, "OLD_SEED")

    os.rename(seed, old_seed)
    cmd = "esl-reformat --mingap stockholm {old_seed} > {seed}".format(old_seed=old_seed, seed=seed)
    subprocess.call(cmd, shell=True)
    if os.path.exists(seed):
        return True

    return False


def parse_arguments():
    parser = argparse.ArgumentParser(description="Precompute families by launching an rfsearch job")
    mutually_exclusive = parser.add_mutually_exclusive_group()
    mutually_exclusive.add_argument("--input", help="A file listing all family directories", type=str)
    mutually_exclusive.add_argument("--csv-input",
                                    help="CSV file with miRNA id, rfam accession number, "
                                         "threshold value of families to update")
    parser.add_argument('--ignore-gap-check', help='Do not remove all gap columns from the alignment', action="store_true",
                        default=False)

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()
    dirs = []
    if args.csv_input:
        rfam_accs = get_rfam_accs(csv_file=args.csv_input)
        for rfam_acc in rfam_accs:
            rfam_family_dir = os.path.join(UPDATE_DIR, rfam_acc)
            dirs.append(rfam_family_dir)
    elif args.input:
        with open(args.input, 'r') as fp:
            dirs = [x.strip() for x in fp]
    if args.ignore_gap_check:
        for family_dir in dirs:
            launch_new_rfsearch(family_dir, cpu=4)
    else:
        for family_dir in dirs:
            if remove_all_gaps(family_dir):
                launch_new_rfsearch(family_dir, cpu=4)
