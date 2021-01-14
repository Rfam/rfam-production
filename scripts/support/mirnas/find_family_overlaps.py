import os
import argparse

from utils import db_utils as db
from scripts.processing import clan_competition as cc

# ----------------------------------------------------------------


def parse_outlist_file(outlist_file):
    """

    :param outlist_file:
    :return:
    """
    outlist_info = {}

    fp = open(outlist_file, 'r')

    for line in fp:
        # if not a comment line
        if line[0] != '#':
            line = line.strip().split()
            unique_accession = "_".join([line[3], line[5], line[6]])
            if unique_accession not in outlist_info:
                outlist_info[unique_accession] = {"evalue": float(line[1]),
                                                  "bit_score": float(line[0]),
                                                  "accession": line[3],
                                                  "start": int(line[5]),
                                                  "end": int(line[6])}

    fp.close()

    return outlist_info


# ----------------------------------------------------------------


def extract_tax_ids_from_species(species_file):
    """

    :param species_file:
    :return:
    """

    tax_ids = {}

    fp = open(species_file, 'r')

    for line in fp:
        # if not a comment line
        if line[0] != '#':
            line = line.strip().split()

            if line[3] not in tax_ids:
                tax_ids[line[3]] = int(line[5])

    fp.close()

    return tax_ids

# ----------------------------------------------------------------


def get_family_location(accession):
    """

    :param family_label:
    :return:
    """

    search_dirs = ["/hps/nobackup/production/xfam/rfam/RELEASES/14.3/miRNA_relabelled/batch1_chunk1_searches",
                   "/hps/nobackup/production/xfam/rfam/RELEASES/14.3/miRNA_relabelled/batch1_chunk2_searches",
                   "/hps/nobackup/production/xfam/rfam/RELEASES/14.3/miRNA_relabelled/batch2/searches"]


    dir_label = ''
    if accession.find("_relabelled") == -1:
        dir_label = accession + "_relabelled"

    for search_dir in search_dirs:
        family_dir_loc = os.path.join(search_dir, dir_label)
        if os.path.exists(family_dir_loc):
            return family_dir_loc

    return None

# ----------------------------------------------------------------

def parse_arguments():
    """

    :return:
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("--family-dir", help="Path to family directory")
    parser.add_argument("--accessions", help="A json file with old/new family mapppings")

    return parser

# ----------------------------------------------------------------

if __name__ == "__main__":

    parser = parse_arguments()
    args = parser.parse_args()

    # load accessions
    fp = open(args.accessions, 'r')
    accessions = json.load(fp)
    fp.close()

    species_file = os.path.join(args.family_dir, "species")
    outlist_file = os.path.join(args.family_dir, "outlist")

    outlist_info = parse_outlist_file(outlist_file)






