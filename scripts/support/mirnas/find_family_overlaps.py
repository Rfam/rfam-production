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


def count_total_num_hits(outlist_hits):

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

    #parser.add_argument("--family-dir", help="Path to family directory")
    parser.add_argument("--accessions", help="A json file with old/new family mapppings")
    parser.add_argument("--add-header", help="Print descriptive header",
                        action="store_true", default=False)
    parser.add_argument("--add-links", help="Creates hyperlinks to available Rfam html content",
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

            outlist_hits = extract_outlist_hits_to_list(outlist_file_loc)

            new_family_unique_accs = list(set(outlist_hits.keys()).difference(set(old_family_full_hits.keys())))
            old_family_unique_accs = list(set(old_family_full_hits.keys()).difference(set(outlist_hits.keys())))

            common_accs = list(set(old_family_full_hits.keys()).intersection(set(outlist_hits.keys())))

            # 3. Find overlaps between the two families
            # only check new family hits
            overlap = 0
            overlap_count = 0

            # count unique counts in intersection
            new_family_count = 0
            old_family_count = 0
            for accession in common_accs:
                old_hits = None
                for region in outlist_hits[accession]:
                    overlap = -1
                    for f_region in old_family_full_hits[accession]:
                        overlap = cc.calc_seq_overlap(region[0], region[1], f_region[0], f_region[1])

                        # this ensures each region in the new family is only checked once
                        # if no overlap found, it iterates over all superfamily hits
                        if overlap >= 0.1:
                           overlap_count += 1

                        # if no overlap detected count old family reqion as unique
                        else:
                            old_family_count += 1

                    # if no overlap count new family region as unique
                    if overlap != -1:
                        new_family_count += 1

            """
            if mirna_id not in family_overlap_counts:
                family_overlap_counts[mirna_id] = {rfam_acc: overlap_count}

            else:
                family_overlap_counts[mirna_id][rfam_acc] = overlap_count
            """

            num_outlist_hits = count_total_num_hits(outlist_hits)

            # calculate old family unique hits
            num_new_family_unique_hits = 0
            for acc in new_family_unique_accs:
                num_new_family_unique_hits += len(outlist_hits[acc])

            # add unique new family regions from intersection
            num_new_family_unique_hits = num_new_family_unique_hits + new_family_count

            # calculate old family unique hits
            num_old_family_unique_hits = 0
            for acc in old_family_unique_accs:
                num_old_family_unique_hits += len(old_family_full_hits[acc])

            # add unique old family regions from intersection
            num_old_family_unique_hits = num_old_family_unique_hits + old_family_count

            total_hits_intersection = 0
            for acc in common_accs:
                total_hits_intersection += len(outlist_hits[acc])
                total_hits_intersection += len(old_family_full_hits[acc])

            # compute family overlap percentage
            overlap_percentage = (float(overlap_count) * 100.0) / float(total_hits_intersection)

            if args.add_header:
                print ("\t".join(["miRBase Id", "Total # new family hits", "# New family unique hits",
                                  "# Overlaps", "Overlaps Percentage",
                                  "# Old family unique hits", "Total # old family hits", "Rfam Acc"]))
                args.add_header = False

            if args.add_links:
                mirbase_hyperlink = "=HYPERLINK(\"https://preview.rfam.org/mirbase/%s_relabelled.html\", \"%s\")"
                mirbase_link = mirbase_hyperlink % (mirna_id, mirna_id)
                mirna_id = mirbase_link

                rfam_hyperlink = "=HYPERLINK(\"https://rfam.org/family/%s\", \"%s\")"
                rfam_link = rfam_hyperlink % (rfam_acc, rfam_acc)
                rfam_acc = rfam_link

            print ("\t".join([mirna_id, str(num_outlist_hits), str(num_new_family_unique_hits),
                              str(overlap_count), str(overlap_percentage), str(num_old_family_unique_hits),
                              total_num_old_family_hits, rfam_acc]))



