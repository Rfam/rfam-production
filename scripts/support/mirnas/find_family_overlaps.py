#!/usr/bin/env python

import os
import json
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

    seen_ga = False

    fp = open(outlist_file, 'r')

    for line in fp:
        # if not a comment line
        if line[0] != '#' and not seen_ga:
            line = line.strip().split()
            unique_accession = "_".join([line[3], line[5], line[6]])
            if unique_accession not in outlist_info:
                outlist_info[unique_accession] = {"evalue": float(line[1]),
                                                  "bit_score": float(line[0]),
                                                  "accession": line[3],
                                                  "start": int(line[5]),
                                                  "end": int(line[6])}

        elif line.find("CURRENT GA THRESHOLD:") != -1:
            seen_ga = True

    fp.close()

    return outlist_info


# ----------------------------------------------------------------


def extract_outlist_hits_to_list(outlist_file, sort=True):
    """

    :param outlist_file:

    :return:
    """

    outlist_hits = {}

    seen_ga = False

    fp = open(outlist_file, 'r')

    for line in fp:
        # if not a comment line
        if line[0] != '#' and not seen_ga:
            line = line.strip().split()

            if line[3] not in outlist_hits:
                outlist_hits[line[3]] = [(int(line[5]), int(line[6]))]
            else:
                outlist_hits[line[3]].append((int(line[5]), int(line[6])))

        elif line.find("CURRENT GA THRESHOLD:") != -1:
            seen_ga = True

    fp.close()

    if sort is True:
        # sorts hits by starting points
        for accession in outlist_hits:
            outlist_hits[accession].sort(key=lambda tup: tup[1])

    return outlist_hits


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


def count_total_outlist_hits(outlist_hits):

    num_hits = 0

    for accession in outlist_hits.keys():
        
        num_hits += len(outlist_hits[accession])

    return num_hits

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

    family_overlap_counts = {}

    for rfam_acc in accessions.keys():
        # fetch family full region hits
        old_family_full_hits = db.fetch_family_full_regions(rfam_acc, sort=True)

        # now work on the miRNA family
        # 1. Detect family dir location
        mirna_id = accessions[rfam_acc]["mirbase_id"]
        family_dir = get_family_location(accessions[rfam_acc]["mirbase_id"])

        # 2. Find outlist path and parse it
        outlist_file_loc = os.path.join(family_dir, "outlist")
        outlist_hits = extract_outlist_hits_to_list(outlist_file_loc)

        # 3. Find overlaps between the two families
        # only check new family hits
        overlap = 0
        overlap_count = 0

        for accession in outlist_hits:
            old_hits = None
            if accession in old_family_full_hits:
                for region in outlist_hits[accession]:
                    for f_region in old_family_full_hits[accession]:
                        overlap = cc.calc_seq_overlap(region[0], region[1], f_region[0], f_region[1])

                        # this ensures each region in the new family is only checked once
                        # if no overlap found, it iterates over all superfamily hits
                        if overlap >= 50:
                            overlap_count += 1
                            break

        family_overlap_counts[mirna_id] = overlap_count
        num_hits = count_total_outlist_hits(outlist_hits)

        # compute family overlap percentage
        overlap_percentage = (family_overlap_counts[mirna_id] * 100) / len(outlist_hits)
        print "\t".join(rfam_acc, mirna_id, str(overlap_percentage))



