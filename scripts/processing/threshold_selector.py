#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
from subprocess import PIPE, Popen
import statistics

# ---------------------------------------------------------------------


def extract_scores_dict_from_outlist_file(scores_file):
    """
    Parses the outlist or species file produced by rfsearch
    and returns all bit_scores above the REVERSED cutoff

    scores_file: The path to rfsearch outlist or species file

    return: A list of all bit scores above REVERSED
    """

    scores = {'SEED': [], 'FULL': [], 'OTHER': []}

    outlist_fp = open(scores_file, 'r')

    for line in outlist_fp:
        if line[0] != '#':
            line = [x for x in line.strip().split(' ') if x!='']
            scores[line[2]].append(float(line[0]))

        else:
            # if we reached REVERSED line, we treat everything as TNs
            # break and return
            if line.find("BEST REVERSED") != -1:
                break

    outlist_fp.close()

    return scores

# ---------------------------------------------------------------------


def extract_bitscores_list_from_scores_file(scores_file):
    """
    Parses the outlist or species file produced by rfsearch
    and returns all bit_scores above the REVERSED cutoff

    scores_file: The path to rfsearch outlist or species file

    return: A list of all bit scores above REVERSED
    """

    scores = []

    outlist_fp = open(scores_file, 'r')

    for line in outlist_fp:
        if line[0] != '#':
            line = [x for x in line.strip().split(' ') if x!='']
            scores.append([float(line[0])])
        else:
            # if we reached REVERSED line, we treat everything as TNs
            # break and return
            if line.find("BEST REVERSED") != -1:
                break

    outlist_fp.close()

    return scores

# ---------------------------------------------------------------------


def is_seed_below_reversed(scores_file):
    """
    Checks if any SEED sequences score below REVERSED

    scores_file: The path to rfsearch outlist or species file

    return: True if SEEDs sequences are found below REVERSED. False
    otherwise
    """

    seen_rev = False

    fp = open(scores_file, 'r')

    for line in fp:
        if line[0] != '#':
            if seen_rev is False:
                continue
            else:
                line = [x for x in line.strip().split('\t') if x!='']
                if line[2] == 'SEED':
                    return True
        elif line.find("REVERSED") != -1:
            seen_rev = True

    fp.close()

    return False

# ---------------------------------------------------------------------


def compute_possible_gathering_thresholds(scores, chunks=4):
    """
    Selects a number of probable gathering thresholds to try
    building the model with.

    param scores:
    param chunks:
    return:
    """

    ga_thresholds = []

    all_scores = scores["SEED"] + scores["FULL"]

    # sorts scores in descending order to match order in outlist
    rev_scores = list(reversed(sorted(all_scores)))

    median = statistics.median(rev_scores)
    min_seed_score = sorted(scores["SEED"])[0]

    index = 0

    if min_seed_score < median:
        ga_thresholds.append(min_seed_score)
        index = rev_scores.index(min_seed_score)
    else:
        if len(rev_scores) % 2 != 0:
            index = rev_scores.index(median)
            ga_thresholds.append(median)
        else:
            # a bit conservative, always chooses the highest value
            index = int(len(rev_scores)/2)
            ga_thresholds.append(rev_scores[index])

    scores_chunk = rev_scores[index:]
    chunks -= 2

    while chunks != 0:
        median = statistics.median(scores_chunk)
        if len(scores_chunk) % 2 != 0:
            index = scores_chunk.index(median)
            ga_thresholds.append(median)
        else:
            # a bit conservative, always chooses the highest value
            index = int(len(scores_chunk)/2)
            ga_thresholds.append(scores_chunk[index])

        scores_chunk = scores_chunk[index:]
        chunks -= 2

    return ga_thresholds

# ---------------------------------------------------------------------


def get_hits_coverage(all_scores, threshold):
    """

    all_scores:

    threshold:

    return:
    """

    # sort scores in descending order
    rev_scores_list = list(reversed(sorted(all_scores)))

    coverage = rev_scores_list.index(threshold)*100/len(all_scores)

    return coverage

# ---------------------------------------------------------------------


def threshold_family_with_rfmake(family_dir, gathering_threshold, full_align=True):
    """
    Calls rfmake.pl to set the gathering threshold of a family

    family_dir: The path to an Rfam family directory
    gathering_threshold: A float value specifying the gathering threshold for the
    family

    return: Full alignment path if it exists or None if not. True upon completion
    and full_align=False
    """

    # change to family directory
    os.chdir(family_dir)

    cmd = "rfmake.pl -t %f"

    if full_align is True:
        cmd = "rfmake.pl -t %f -a"

    subprocess.call(cmd % gathering_threshold, shell=True)

    if full_align is True:
        full_align_path = os.path.join(family_dir, "align")
        if os.path.exists(full_align_path):
            return full_align_path
        else:
            return None

    return True

# ---------------------------------------------------------------------


def generate_family_ss_with_rscape(family_dir, file_type='SEED'):
    """

    family_dir:
    file_type: The type of the alignment SEED/FULL

    return:
    """

    alignment_path = os.path.join(family_dir, file_type)
    
    if file_type == "FULL":
        alignment_path = os.path.join(family_dir, "align")

    outdir = os.path.join(family_dir, "rscape-" + file_type.lower())

    # create outdir if it does not exist
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    cmd = "R-scape --outdir %s -s --fold %s" % (outdir, alignment_path)

    process = Popen(["R-scape", "--outdir", outdir,
                     "-s", "--fold", alignment_path], stdout=PIPE)

    rscape_out = process.communicate()[0].decode('Utf-8').split('\n')

    covarying_bp = rscape_out[-3].split(' ')[-1]
    
    return int(covarying_bp)

# ---------------------------------------------------------------------


def parse_arguments():
    """
    Basic argument parsing using Python's argparse

    return: A valid argpase parser object
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("--family-dir", help='The path to an Rfam family directory',
                        action='store')

    parser.add_argument("--multi", help='Specifies a directory with multiple Rfam family directories',
                        action='store_true')

    return parser

# ---------------------------------------------------------------------


if __name__ == '__main__':

    parser = parse_arguments()

    args = parser.parse_args()

    max_covariation = 0

    outlist = os.path.join(args.family_dir, "outlist")

    scores = extract_scores_dict_from_outlist_file(outlist)

    threshold_list = compute_possible_gathering_thresholds(scores, chunks=6)

    best_threshold = (-1, -1)
    
    for bit_score in threshold_list:
        full_path = threshold_family_with_rfmake(args.family_dir, bit_score, full_align=True)

        if not os.path.exists(full_path):
            sys.exit("Error generating FULL alignment for family %s" % os.path.basename(args.family_dir))

        covarying_bp = generate_family_ss_with_rscape(args.family_dir, file_type='FULL')

        if covarying_bp > max_covariation:
            max_covariation = covarying_bp
            best_threshold = (bit_score, covarying_bp)

    full_path = threshold_family_with_rfmake(args.family_dir, best_threshold[0], full_align=True)

