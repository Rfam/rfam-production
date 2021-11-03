import json
import os
import argparse
import subprocess
from scripts.support.mirnas.update_mirnas_helpers import MEMORY, CPU, LSF_GROUP, get_rfam_accs
from scripts.support.mirnas.config import UPDATE_DIR, SEARCH_DIRS


def add_ref_sequentially(reference, thresholds_file=None, rfam_accessions=None):
    """
    Call add_ref.pl sequentially
    :param reference: PubMed reference ID for the DESC, by default the latest MiRBase paper 30423142
    :param rfam_accessions: list of Rfam accession numbers
    :param thresholds_file: JSON file with miRNA IDs and the corresponding threshold
    """
    cmd = ("add_ref.pl {0}".format(reference))

    def call_addref(fam_dir):
        if os.path.exists(fam_dir):
            os.chdir(fam_dir)
            subprocess.call(cmd, shell=True)

    if rfam_accessions:
        for family in rfam_accessions:
            family_dir = os.path.join(UPDATE_DIR, family)
            call_addref(family_dir)

    elif thresholds_file:
        with open(thresholds_file, 'r') as fp:
            thresholds = json.load(fp)
        for family in thresholds.keys():
            for searchdir in SEARCH_DIRS:
                if family.find("relabelled") == -1:
                    family_dir = os.path.join(searchdir, family + "_relabelled")
                else:
                    family_dir = os.path.join(searchdir, family)
                call_addref(family_dir)


def call_add_ref_cmd(fam_dir, ref):
    """
    Submit the add_ref.pl jpb
    :param fam_dir: family directory, from which to call the command
    :param ref: PubMed reference ID for the DESC, by default the latest MiRBase paper 30423142
    """
    cmd = ("bsub -M {mem} -o {out_file} -e {err_file} -n {cpu} -g {lsf_group} -q production-rh74 "
           " \"cd {family_dir} && add_ref.pl {ref}\"")
    lsf_err_file = os.path.join(fam_dir, "auto_add_ref.err")
    lsf_out_file = os.path.join(fam_dir, "auto_add_ref.out")
    subprocess.call(
        cmd.format(
            mem=MEMORY, out_file=lsf_out_file, err_file=lsf_err_file, cpu=CPU, lsf_group=LSF_GROUP, family_dir=fam_dir,
            ref=ref), shell=True)


def auto_add_ref(reference, rfam_accessions=None, thresholds_file=None):
    """
    Call add_ref.pl
    :param reference: PubMed reference ID for the DESC, by default the latest MiRBase paper 30423142
    :param rfam_accessions: list of Rfam accession numbers
    :param thresholds_file: JSON file with miRNA IDs and the corresponding threshold
    """
    if rfam_accessions:
        for family in rfam_accessions:
            family_dir = os.path.join(UPDATE_DIR, family)
            if os.path.exists(family_dir):
                call_add_ref_cmd(family_dir, reference)
            else:
                continue
    elif thresholds_file:
        with open(thresholds_file, 'r') as fp:
            thresholds = json.load(fp)

        for family in thresholds.keys():
            for searchdir in SEARCH_DIRS:
                if family.find("relabelled") == -1:
                    family_dir = os.path.join(searchdir, family + "_relabelled")
                else:
                    family_dir = os.path.join(searchdir, family)

                if os.path.exists(family_dir):
                    call_add_ref_cmd(family_dir, reference)
                else:
                    continue


def parse_arguments():
    parser = argparse.ArgumentParser()
    file_input = parser.add_mutually_exclusive_group()
    file_input.add_argument("--csv-input",
                            help="CSV file with miRNA id, rfam accession number, threshold value of families to update")
    file_input.add_argument("--mirna-list", help="A json file with all miRNA candidates",
                            action="store")
    parser.add_argument("--ref", help="A string indicating the PubMed id to use for reference",
                        action="store", default="30423142")
    parser.add_argument("--sequential", help="Modify DESC files sequentially",
                        action="store_true", default=False)

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()
    if args.csv_input:
        rfam_accs = get_rfam_accs(args.csv_input)
    else:
        rfam_accs = None
    if args.mirna_list:
        mirna_list = args.mirna_list
    else:
        mirna_list = None
    if args.sequential is False:
        auto_add_ref(reference=args.ref, rfam_accessions=rfam_accs, thresholds_file=mirna_list)
    else:
        add_ref_sequentially(reference=args.ref, rfam_accessions=rfam_accs, thresholds_file=mirna_list)
