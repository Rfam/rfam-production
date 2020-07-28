import os
import sys
import argparse
import subprocess
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
            index = (len(rev_scores)/2)-1
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
            index = (len(rev_scores) / 2) - 1
            ga_thresholds.append(rev_scores[index])

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

    outlist = "../../data/outlist"

    scores = extract_scores_dict_from_outlist_file(outlist)
    all_scores = scores['SEED'] + scores['FULL']

    rev_scores_list = list(reversed(sorted(all_scores)))
    print (rev_scores_list)

    print (compute_possible_gathering_thresholds(scores, chunks=6))

    """
    search_dir = sys.argv[1]

    family_dirs = os.listdir(search_dir)

    for dir in family_dirs:
        dir_loc = os.path.join(search_dir, dir)
        outlist_loc = os.path.join(dir_loc, )

    if is_seed_below_reversed(outlist) is True:
        sys.exit("ERROR: SEED sequences below reversed!!")

    scores = extract_scores_dict_from_outlist_file(outlist)

    all_scores = []
    all_scores = scores['SEED'] + scores['FULL']

    #print ("SEEDs: ", scores['SEED'])
    #print ("FULLs: ", scores['FULL'])
    #print ("all_scores: ", sorted(all_scores))
    min_seed = sorted(scores['SEED'])[0]
    median = statistics.median(all_scores)
    num_scores = len(all_scores)


    print ("median: ", statistics.median(all_scores))
    print ("avg: ", statistics.mean(all_scores))
    print ("stdev: ", statistics.stdev(all_scores))
    print ("#scores: ", len(all_scores))
    print ("index: ", sorted(all_scores).index(52.5))

    rev_scores = list(reversed(sorted(all_scores)))
    print ("rev_scores: ", rev_scores)
    print ("threshold coverage: ", rev_scores.index(52.5)*100/len(all_scores))

    """

