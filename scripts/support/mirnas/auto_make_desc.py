import os
import sys
import json
import argparse
import subprocess

searchdirs = ["/hps/nobackup/production/xfam/rfam/RELEASES/14.3/miRNA_relabelled/batch1_chunk1_searches",
              "/hps/nobackup/production/xfam/rfam/RELEASES/14.3/miRNA_relabelled/batch1_chunk2_searches",
              "/hps/nobackup/production/xfam/rfam/RELEASES/14.3/miRNA_relabelled/batch2/searches"]

MEMORY = 10000
CPU = 4
LSF_GROUP = "/family_srch"

env_path = "/nfs/production/xfam/users/rfamprod/code/env2/bin/activate"
desc_generator_path = "/nfs/production/xfam/users/rfamprod/code/rfam-production/scripts/preprocessing/desc_generator.py"


# --------------------------------------------------------------------------------------------------

def auto_desc_make(family_accession, wiki_links_file):
    cmd = ("bsub -M %s -o %s -e %s -n %s -g %s -q production-rh74 "
           "-J %s \"source %s && python %s --input %s --wiki-links %s\"")

    for searchdir in searchdirs:
        family_dir = ""
        if family_accession.find("relabelled") == -1:
            family_dir = os.path.join(searchdir, family_accession + "_relabelled")
        else:
            family_dir = os.path.join(searchdir, family_accession)

        if os.path.exists(family_dir):
            lsf_err_file = os.path.join(family_dir, "auto_desc_make.err")
            lsf_out_file = os.path.join(family_dir, "auto_desc_make.out")
            job_name = family_accession.partition('_')[0]

            # print cmd % (MEMORY, lsf_out_file, lsf_err_file, CPU, LSF_GROUP, job_name, family_dir, str(thresholds[family]))
            subprocess.call(cmd % (MEMORY, lsf_out_file, lsf_err_file, CPU, LSF_GROUP, job_name, env_path,
                                   desc_generator_path, family_dir, wiki_links_file), shell=True)

        else:
            continue

    return family_dir


# --------------------------------------------------------------------------------------------------


def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument("--mirna-ids", help="A .json file containing all miRNAs to generate DESC files for",
                        action="store")
    parser.add_argument("--wiki-links", help="A file listing all possible miRNA wiki links",
                        action="store")
    parser.add_argument("--no-ref", help="Rewrites the DESC omitting any reference lines",
                        action="store_true", default=False)
    return parser


# --------------------------------------------------------------------------------------------------


if __name__ == '__main__':

    parser = parse_arguments()
    args = parser.parse_args()

    json_file = args.mirna_ids
    wiki_links = args.wiki_links

    fp = open(json_file, 'r')
    accessions = json.load(fp)
    fp.close()

    for miRNA_family in accessions.keys():
        # print (miRNA_family)
        if miRNA_family[0:4] == "MIPF":
            family_dir = auto_desc_make(miRNA_family, wiki_links)
            # print (family_dir)
            # sys.exit()
