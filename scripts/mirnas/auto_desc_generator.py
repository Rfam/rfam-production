import os
import argparse
import subprocess

from scripts.mirnas.mirna_config import NEW_DIR, MEMORY, CPU, LSF_GROUP, ENV_PATH, DESC_GEN_PATH
from scripts.mirnas.update_mirnas_helpers import get_mirna_ids


def auto_desc_make(mirna_id, wiki_links_file=None):
    """
    Launch the jobs to create the new DESC file

    :param mirna_id: miRBase family ID
    :param wiki_links_file: A file listing all possible miRNA wiki links
    """
    family_dir = os.path.join(NEW_DIR, mirna_id)
    if os.path.exists(family_dir):
        lsf_err_file = os.path.join(family_dir, "auto_desc_make.err")
        lsf_out_file = os.path.join(family_dir, "auto_desc_make.out")

        if wiki_links_file:
            cmd = ("bsub -M {mem} -o {out_file} -e {err_file} -n {cpu} -g {lsf_group} "
                   "-J {job_name} \"source {env_path} && python {desc_gen} --input {family_dir} --wiki-links {wiki}\"")
            subprocess.call(
                cmd.format(
                    mem=MEMORY, out_file=lsf_out_file, err_file=lsf_err_file, cpu=CPU, lsf_group=LSF_GROUP,
                    job_name=mirna_id, env_path=ENV_PATH, desc_gen=DESC_GEN_PATH, family_dir=family_dir,
                    wiki=wiki_links_file), shell=True)
        else:
            cmd = ("bsub -M {mem} -o {out_file} -e {err_file} -n {cpu} -g {lsf_group} "
                   "-J {job_name} \"source {env_path} && python {desc_gen} --input {family_dir}\"")
            subprocess.call(
                cmd.format(
                    mem=MEMORY, out_file=lsf_out_file, err_file=lsf_err_file, cpu=CPU, lsf_group=LSF_GROUP,
                    job_name=mirna_id, env_path=ENV_PATH, desc_gen=DESC_GEN_PATH, family_dir=family_dir,), shell=True)

    else:
        print('error')


def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument("--input", help="A .json file mirna : threshold pairs",
                        action="store")
    parser.add_argument("--wiki-links", help="A file listing all possible miRNA wiki links",
                        action="store")
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()

    wiki_links = args.wiki_links
    mirna_ids = get_mirna_ids(args.input)
    for mirna in mirna_ids:
        print(mirna)
        auto_desc_make(mirna, wiki_links)
