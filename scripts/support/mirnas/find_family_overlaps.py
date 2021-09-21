#!/usr/bin/env python

import os
import json
import argparse
import re

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


def extract_outlist_hits_to_dict(outlist_file, skip_seed=True, sort=True):
    """

    :param outlist_file:

    :return:
    """
    outlist_hits = {}
    seen_ga = False

    fp = open(outlist_file, 'r')

    for line in fp:
        # skip SEED sequences if skip_seed option is enabled
        if skip_seed is True and 'SEED' in line:
            continue
        # if not a comment line
        if line[0] != '#' and not seen_ga:
            line = line.strip().split()

            if line[3] not in outlist_hits:
                outlist_hits[line[3]] = [(int(line[5]), int(line[6]))]
            else:
                outlist_hits[line[3]].append((int(line[5]), int(line[6])))

        elif 'CURRENT GA THRESHOLD' in line:
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

def compare_ids(mirbase_id, rfam_id):
    """
    Compare ids like MIPF0000024__mir-103 and mir-103
    """
    parts = mirbase_id.split('__')
    if parts[-1].lower() == rfam_id.lower():
        return 'Yes'
    else:
        return 'No'

# ----------------------------------------------------------------


def parse_overlaps(overlaps_file):
    """
    External overlap [SS] of CM000298.1/96020656-96020744 with RF00998:CM000298.1/96020650-96020744 by fullOL
    External overlap [SS] of CM000316.3/123249182-123249270 with RF00998:CM000316.3/123249176-123249270 by fullOL
    """
    overlap_accs = set()
    with open(overlaps_file, 'r') as f:
        for line in f:
            if 'External overlap' not in line:
                continue
            match = re.search(r'with (RF\d{5})\:', line)
            if match:
                overlap_accs.add(match.group(1))
    return overlap_accs
# ----------------------------------------------------------------

def get_desc_ga(filename):
    """
    Get the gathering threshold (GA) from the DESC file.
    """
    with open(filename, 'r') as f:
        for line in f:
            if not line.startswith('GA'):
                continue
            parts = re.split(r'\s+', line.strip())
            return float(parts[-1])


def main(args):
    # load accessions
    accessions = None
    with open(args.accessions, 'r') as fp:
        accessions = json.load(fp)

    row_id = 3 # first two rows are taken by the header
    # iterate over all miRBase ids
    for mirna_id in accessions.keys():

        # 1. Detect family dir location
        family_dir = get_family_location(mirna_id)
        if 'MIPF0000482__mir' in mirna_id:
            continue
        # if 'MIPF0000858__mir' not in mirna_id:
        #     continue

        # check if the threshold has been set correctly
        desc_ga = get_desc_ga(os.path.join(family_dir, 'DESC'))
        curated_ga = float(accessions[mirna_id]['threshold'][0])
        if desc_ga != curated_ga:
            cmd = 'cd {} && rfmake.pl -t {} && cd -'.format(family_dir, curated_ga)
            print(cmd)
            os.system(cmd)

        # calculate overlaps
        overlaps = os.path.join(family_dir, 'overlap')
        if not os.path.exists(overlaps):
            parent_dir = os.path.dirname(family_dir)
            cmd = 'cd {} && rqc-overlap.pl {} && cd -'.format(parent_dir, os.path.basename(family_dir))
            os.system(cmd)
        rfam_accs = parse_overlaps(overlaps)

        # iterate over all rfam_accs overlapping with new miRBase candidates
        for rfam_acc in rfam_accs:

            # fetch family full region hits
            old_family_full_hits = db.fetch_family_full_regions(rfam_acc, sort=True)
            total_num_old_family_hits = count_total_num_hits(old_family_full_hits)

            # 2. Find outlist path and parse it
            outlist_file_loc = os.path.join(family_dir, "outlist")

            # skip family if the outlist does not exist - possible search job crash
            if not os.path.exists(outlist_file_loc):
                print('Warning: an outlist file not found {}'.format(outlist_file_loc))
                continue

            # extract new family hits from outlist file
            outlist_hits = extract_outlist_hits_to_dict(outlist_file_loc, skip_seed=True)
            # total number of hits extracted from the outlist filen
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
            old_family_overlapped = set()
            for accession in common_accs:
                for region in outlist_hits[accession]:
                    found_overlap = False
                    for f_region in old_family_full_hits[accession]:
                        if cc.get_strand(region[0], region[1]) != cc.get_strand(f_region[0], f_region[1]):
                            continue
                        overlap = cc.calc_seq_overlap(region[0], region[1], f_region[0], f_region[1])
                        if overlap > 0:
                            overlap_count += 1
                            found_overlap = True
                            old_family_overlapped.add('{}_{}_{}'.format(accession, f_region[0], f_region[1]))
                    # if no overlap count new family region as unique
                    if not found_overlap:
                        new_unique_intersect += 1

            # check if there are any old hits that did not overlap anything
            for accession in common_accs:
                for f_region in old_family_full_hits[accession]:
                    key = '{}_{}_{}'.format(accession, f_region[0], f_region[1])
                    if key not in old_family_overlapped:
                        old_unique_intersect += 1

            # Now add the overlap counts
            # add unique new family regions from intersection
            num_new_family_unique_hits = num_new_family_unique_hits + new_unique_intersect

            # add unique old family regions from intersection
            num_old_family_unique_hits = num_old_family_unique_hits + old_unique_intersect

            rfam_acc_metadata = db.fetch_family_metadata(rfam_acc)

            if args.add_links:
                mirbase_hyperlink = "=HYPERLINK(\"https://preview.rfam.org/mirbase/%s_relabelled.html\", \"%s\")"
                mirbase_link = mirbase_hyperlink % (mirna_id, mirna_id)

                rfam_hyperlink = '=HYPERLINK("https://rfam.org/family/{0}", "{0}")'
                rfam_link = rfam_hyperlink.format(rfam_acc)
            else:
                mirbase_link = mirna_id
                rfam_link = rfam_acc

            # tax id variable declaration
            new_family_taxids = None
            old_family_taxids = None

            if args.check_taxids:
                species_file = os.path.join(family_dir, "species")
                # find which tax id is missing from essential species
                new_family_taxids = [str(x) for x in list(set(ESSENTIAL_TAXIDS).intersection(set(extract_tax_ids_from_species_file(species_file).values())))]
                old_family_taxids = [str(x) for x in list(set(ESSENTIAL_TAXIDS).intersection(set(db.fetch_family_tax_ids(rfam_acc))))]

                new_family_taxids_str = "N/A"
                old_family_taxids_str = "N/A"

                if len(new_family_taxids) == 1:
                    new_family_taxids_str = new_family_taxids[0]
                elif len(new_family_taxids) == 2:
                    new_family_taxids_str = ','.join(new_family_taxids)

                if len(old_family_taxids) == 1:
                    old_family_taxids_str = old_family_taxids[0]
                elif len(old_family_taxids) == 2:
                    old_family_taxids_str = ','.join(old_family_taxids)

            print ("\t".join([
                mirbase_link,
                str(total_num_outlist_hits),
                str(num_new_family_unique_hits),
                str(overlap_count),
                str(num_old_family_unique_hits),
                str(total_num_old_family_hits),
                rfam_link,
                rfam_acc_metadata['rfam_id'],
                '{:.1f}'.format((float(overlap_count) * 100)/total_num_outlist_hits),
                '{:.1f}'.format((float(overlap_count) * 100)/total_num_old_family_hits),
                compare_ids(mirna_id, rfam_acc_metadata['rfam_id']),
                new_family_taxids_str,
                old_family_taxids_str,
                'if(AND(J{0}>=85,K{0}="Yes"), "Update", "Review")'.format(row_id),
            ]))
            row_id += 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--accessions", help="A json file with old/new family mappings")
    parser.add_argument("--add-links", help="Creates hyperlinks to available Rfam html content",
                        action="store_true", default=False)
    parser.add_argument("--check-taxids", help="Check essential taxonomy ids exist",
                        action="store_true", default=False)
    main(parser.parse_args())
