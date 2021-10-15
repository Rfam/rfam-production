import os
import subprocess
import argparse

from scripts.support.mirnas.update_mirnas_helpers import MEMORY, LSF_GROUP, get_rfam_accs
from scripts.support.mirnas.config import UPDATE_DIR


def launch_new_rfsearch(family_dir, cpu):
    """
    Launches a new LSF job

    family_dir: The location of a valid family directory (directory with Rfam ID)
    cpus: Number of CPUs to use per thread
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
    seed = os.path.join(family_dir, "SEED")
    old_seed = os.path.join(family_dir, "OLD_SEED")

    # rename SEED to OLD_SEED
    os.rename(seed, old_seed)
    # regenerate SEED file removing all gaps
    cmd = "esl-reformat --mingap stockholm %s > %s" % (old_seed, seed)

    subprocess.call(cmd, shell=True)

    if os.path.exists(seed):
        return True

    return False


def parse_arguments():
    parser = argparse.ArgumentParser(description='Script to precompute families')
    required_arguments = parser.add_argument_group("required arguments")
    required_arguments.add_argument("--input", help="CSV file", type=str)
    parser.add_argument('--no-gap', help='Remove all gap columns from the alignment', action="store_true",
                        default=False)

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()
    rfam_accs = get_rfam_accs(csv_file=args.input)
    for rfam_acc in rfam_accs:
        # commenting out for now as I don't know if it's needed
        # if args.no_gap:
        #     check = remove_all_gaps(rfam_family_dir)
        #     if check is True:
        rfam_family_dir = os.path.join(UPDATE_DIR, rfam_acc)
        launch_new_rfsearch(rfam_family_dir, cpu=4)
