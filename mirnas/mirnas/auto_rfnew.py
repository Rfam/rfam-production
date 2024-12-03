import argparse
import os
import subprocess

from scripts.mirnas.mirna_config import CPU, LSF_GROUP, MEMORY, NEW_DIR
from scripts.mirnas.update_mirnas_helpers import get_mirna_ids


def auto_rfnew(mirnas):
    """
    Run rfnew.pl to check-in the given families
    :param mirnas: list of miRNA IDs to check-in
    """
    for mirna in mirnas:
        family_dir = os.path.join(NEW_DIR, mirna)
        lsf_err_file = os.path.join(family_dir, "auto_rfnew.err")
        lsf_out_file = os.path.join(family_dir, "auto_rfnew.out")
        cmd = (
            "bsub -M {mem} -o {out_file} -e {err_file} -n {cpu} -g {lsf_group} "
            "-J {job_name} \"cd {new_dir} && rfnew.pl -m 'Adding new miRNA family' {mirna_id}\""
        )
        subprocess.call(
            cmd.format(
                mem=MEMORY,
                out_file=lsf_out_file,
                err_file=lsf_err_file,
                cpu=CPU,
                lsf_group=LSF_GROUP,
                job_name=mirna,
                new_dir=NEW_DIR,
                mirna_id=mirna,
            ),
            shell=True,
        )


def parse_arguments():
    parser = argparse.ArgumentParser(description="run rfnew.pl")
    parser.add_argument(
        "--input",
        help="TSV file with miRNA ID, and threshold value of families to update",
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()
    mirna_ids = get_mirna_ids(args.input)
    auto_rfnew(mirna_ids)
