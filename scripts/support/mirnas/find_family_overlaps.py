#!/usr/bin/env python

import os
import json
import argparse

from utils import db_utils as db
from scripts.processing import clan_competition as cc

# ----------------------------------------------------------------

# human and mouse tax ids
ESSENTIAL_TAXIDS = [9606, 10090]

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


def extract_outlist_hits_to_dict(outlist_file, sort=True):
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


def extract_tax_ids_from_species_file(species_file):
    """
    Parses family's species file and extracts all distinct tax ids

    :param species_file: The path to a family's species file

    :return: A dictionary of
    """

    tax_ids = {}
    seen_ga = False

    fp = open(species_file, 'r')

    for line in fp:
        # if not a comment line
        if line[0] != '#' and not seen_ga:
            line = line.strip().split()
            
            if line[3] not in tax_ids:
                if line[5] != '-':
                    tax_ids[line[3]] = int(line[5])
        
        elif line.find("CURRENT GA THRESHOLD:") != -1:
            seen_ga = True

    fp.close()

    return tax_ids

# ----------------------------------------------------------------


def get_family_location(family_label):
    """
    Finds exact location of a miRNA familly based on the miRBase
    accession and the destination directory

    :param family_label: Family directory label

    :return: A directory location if it exists, None otherwise
    """

    search_dirs = ["/hps/nobackup/production/xfam/rfam/RELEASES/14.3/miRNA_relabelled/batch1_chunk1_searches",
                   "/hps/nobackup/production/xfam/rfam/RELEASES/14.3/miRNA_relabelled/batch1_chunk2_searches",
                   "/hps/nobackup/production/xfam/rfam/RELEASES/14.3/miRNA_relabelled/batch2/searches"]

    dir_label = ''
    if family_label.find("_relabelled") == -1:
        dir_label = family_label + "_relabelled"

    for search_dir in search_dirs:
        family_dir_loc = os.path.join(search_dir, dir_label)
        if os.path.exists(family_dir_loc):
            return family_dir_loc

    return None

# ----------------------------------------------------------------


def count_total_num_hits(outlist_hits):
    """
    Counts total number of family hits

    :param outlist_hits: A dictionary in the form of {rfamseq_acc: [(s1,e1),...]

    :return: Total number of hits found in the dictionary
    """

    total_num_hits = 0

    for accession in outlist_hits.keys():
        total_num_hits += len(outlist_hits[accession])

    return total_num_hits

# ----------------------------------------------------------------


def parse_arguments():
    """
    Basic argument parsing using Argparse

    :return: An argparse parser object
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("--accessions", help="A json file with old/new family mapppings")
    parser.add_argument("--add-header", help="Print descriptive header",
                        action="store_true", default=False)
    parser.add_argument("--add-links", help="Creates hyperlinks to available Rfam html content",
                        action="store_true", default=False)
    parser.add_argument("--check-taxids", help="Check essential taxonomy ids exist",
                        action="store_true", default=False)

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

    # iterate over all miRBase ids
    for mirna_id in accessions.keys():

        rfam_accs = accessions[mirna_id]["rfam_acc"]

        # iterate over all rfam_accs overlapping with new miRBase candidates
        for rfam_acc in rfam_accs:

            # skip iteration if not a valid Rfam family accession
            if rfam_acc[0:2] != "RF":
                continue

            # fetch family full region hits
            old_family_full_hits = db.fetch_family_full_regions(rfam_acc, sort=True)
            total_num_old_family_hits = count_total_num_hits(old_family_full_hits)

            # now work on the miRNA family
            # 1. Detect family dir location
            family_dir = get_family_location(mirna_id)

            # 2. Find outlist path and parse it
            outlist_file_loc = os.path.join(family_dir, "outlist")

            # skip family if the outlist does not exist - possible search job crash
            if not os.path.exists(outlist_file_loc):
                continue

            # extract new family hits from outlist file
            outlist_hits = extract_outlist_hits_to_dict(outlist_file_loc)
            # total number of hits extracted from the outlist file
            total_num_outlist_hits = count_total_num_hits(outlist_hits)

            # calculate new family unique hits
            # 1. find new family unique accessions
            new_family_unique_accs = list(set(outlist_hits.keys()).difference(set(old_family_full_hits.keys())))

            # 2. count number of hits per unique accession belonging to family
            num_new_family_unique_hits = 0
            for acc in new_family_unique_accs:
                num_new_family_unique_hits += len(outlist_hits[acc])

            #######

            # calculate old family unique hits
            # 1. find old family unique accessions
            old_family_unique_accs = list(set(old_family_full_hits.keys()).difference(set(outlist_hits.keys())))

            # 2. count number of hits per unique accession belonging to family
            num_old_family_unique_hits = 0
            for acc in old_family_unique_accs:
                num_old_family_unique_hits += len(old_family_full_hits[acc])

            # Now work on finding overlaps between the two families
            # 1. find any common accessions between the families
            common_accs = list(set(old_family_full_hits.keys()).intersection(set(outlist_hits.keys())))

            # 2. Find overlaps between hits of the two families
            # only checks new family hits
            overlap = 0
            overlap_count = 0

            # count unique counts in intersection
            new_unique_intersect = 0
            old_unique_intersect = 0
            for accession in common_accs:
                old_hits = None
                new_unique_hit = False
                for region in outlist_hits[accession]:
                    overlap = -1
                    found_overlap = False
                    for f_region in old_family_full_hits[accession]:
                        overlap = cc.calc_seq_overlap(region[0], region[1], f_region[0], f_region[1])

                        # this ensures each region in the new family is only checked once
                        # if no overlap found, it iterates over all superfamily hits
                        if overlap > 0:
                           overlap_count += 1
                           found_overlap = True

                        # if no overlap detected count old family reqion as unique
                        else:
                            old_unique_intersect += 1
                            found_overlap = True

                    # if no overlap count new family region as unique
                    if not found_overlap:
                        new_unique_intersect += 1

            # Now add the overlap counts
            # add unique new family regions from intersection
            num_new_family_unique_hits = num_new_family_unique_hits + new_unique_intersect

            # add unique old family regions from intersection
            num_old_family_unique_hits = num_old_family_unique_hits + old_unique_intersect

            if args.add_header:
                print ("\t".join(["miRBase Id", "Total # new family hits", "# New family unique hits",
                                  "# Overlaps", "# Old family unique hits", "Total # old family hits",
                                  "Rfam Acc"]))
                args.add_header = False

            if args.add_links:
                mirbase_hyperlink = "=HYPERLINK(\"https://preview.rfam.org/mirbase/%s_relabelled.html\", \"%s\")"
                mirbase_link = mirbase_hyperlink % (mirna_id, mirna_id)
                mirna_id = mirbase_link

                rfam_hyperlink = "=HYPERLINK(\"https://rfam.org/family/%s\", \"%s\")"
                rfam_link = rfam_hyperlink % (rfam_acc, rfam_acc)
                rfam_acc = rfam_link

	    if args.check_taxids:
		pass

            print ("\t".join([mirna_id, str(total_num_outlist_hits), str(num_new_family_unique_hits),
                              str(overlap_count), str(num_old_family_unique_hits),
                              str(total_num_old_family_hits), rfam_acc]))
