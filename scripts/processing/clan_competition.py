#!/usr/bin/python
"""
Copyright [2009-2016] EMBL-European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
     http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

"""
Description: Clan competition script

Notes: clan files are generated using export script clan_file_generator.py
       and sorted on rfamseq_acc (col2) using linux sort command as:
       sort -k2 -t $'\t\' clan_file.txt > clan_file_sorted.txt
"""

# ---------------------------------IMPORTS-------------------------------------

import sys
import os
import string
import logging
import timeit
from utils import db_utils

# -----------------------------------------------------------------------------

RFAM_ACC = 0  # full region rfam_acc
SEQ_ACC = 1  # full region rfamseq_acc
START = 2  # full region seq_start
END = 3  # full region seq_end
EVAL = 5  # full region evalue score
TRUNC = 8  # full region truncated

OVERLAP = 0.5  # overlap cutoff
COMP_OVL = 1.0  # complete overlap
NO_OVL = 0.0    # no overlap

# -----------------------------------------------------------------------------


def get_strand(start, end):
    """
    Checks the start and end coordinates of a sequence and returns -1 if the
    sequence comes from the 3' strand and 1 if it comes from the 5' strand

    start: An integer indicating sequence start location
    end: An integer indicating sequence end location
    """

    # -
    if start > end:
        return -1
    # +
    elif start <= end:
        return 1

    return 0

# -----------------------------------------------------------------------------


def calc_seq_overlap(s1, e1, s2, e2):
    """
    Calculate sequence overlaps

    s1: SEQ1 start coordinate
    e1: SEQ1 end coordinate
    s2: SEQ2 start coordinate
    e2: SEQ2 end coordinate
    """

    len1 = abs(e1 - s1)
    len2 = abs(e2 - s2)

    overlap = None

    # get strand
    strand = get_strand(s1, e1)

    # full overlap
    if s1 == s2 and len1 == len2:
        return float(len1) / float(len2)

    # check5'
    elif strand == 1:
        overlap = cal_overlap_pos_strand(s1, e1, s2, e2)

    # check 3'
    elif strand == -1:
        overlap = cal_overlap_neg_strand(s1, e1, s2, e2)

    # will return None in a case that we didn't capture
    return overlap

# -----------------------------------------------------------------------------


def cal_overlap_pos_strand(s1, e1, s2, e2):
    """
    Calculates the region overlap between two regions on the 5' strand and
    returns the degree of overlap

    s1: Seq1 start coordinates
    e1: Seq1 end coordinates
    s2: Seq2 start coordinates
    e2: Seq2 end coordinates
    """

    overlap = None

    len1 = abs(e1 - s1)
    len2 = abs(e2 - s2)

    min_len = min(len1, len2)

    # seq2 within seq1
    if s1 < s2 and e2 < e1:
        overlap = COMP_OVL

    # partial overlap, seq1 before seq2
    elif s1 <= s2 and s2 < e1 and e1 <= e2:
        overlap = float(e1 - s2 + 1) / float(min_len)

    # no overlap, seq1 before seq2
    elif s1 < s2 and e1 <= s2:
        overlap = NO_OVL

    # seq1 within seq2 region
    elif s2 < s1 and e1 < e2:
        overlap = COMP_OVL

    # no overlap, seq2 before seq1
    elif s2 < s1 and e2 <= s1:
        overlap = NO_OVL

    # partial overlap, seq2 before seq1
    elif s2 <= s1 and s1 < e2 and e2 <= e1:
        overlap = float(e2 - s1 + 1) / float(min_len)

    return overlap

# -----------------------------------------------------------------------------


def cal_overlap_neg_strand(s1, e1, s2, e2):
    """
    Calculates the region overlap between two regions on the 3' strand and
    returns the degree of overlap

    s1: Seq1 start coordinates
    e1: Seq1 end coordinates
    s2: Seq2 start coordinates
    e2: Seq2 end coordinates

    """

    overlap = None

    len1 = abs(e1 - s1)
    len2 = abs(e2 - s2)

    min_len = min(len1, len2)

    # seq2 within seq1 region - this may match the partial overlap case
    if s1 > s2 and e1 < e2:
        overlap = COMP_OVL

    # no overlap, seq1 before seq2
    elif s1 > s2 and e1 >= s2:
        overlap = NO_OVL

    # partial overlap, seq1 before seq2
    elif s1 >= s2 and s2 > e1 and e1 >= e2:
        overlap = float(s2 - e1 + 1) / float(min_len)

    # seq1 within seq2 region
    elif s2 > s1 and e1 > e2:
        overlap = COMP_OVL

    # no overlap, seq2 before seq1
    elif s2 > s1 and e2 >= s1:
        overlap = NO_OVL

    # partial overlap, seq2 before seq1
    elif s2 >= s1 and s1 > e2 and e2 >= e1:
        overlap = float(s1 - e2 + 1) / float(min_len)

    return overlap

# -----------------------------------------------------------------------------


def compete_seq_regions(regions, log):
    """
    regions: A list of duplicate regions for seq_acc
    log: log file pointer for tracking regions we haven't captured
    """

    index = 0

    non_sig_regs = []

    while index <= len(regions) - 2:
        reg1 = regions[index]
        comp_regs = regions[index + 1:]

        for reg2 in comp_regs:

            strand1 = get_strand(int(reg1[START]), int(reg1[END]))
            strand2 = get_strand(int(reg2[START]), int(reg2[END]))

            # check if the sequences come from the same strand
            if strand1 == strand2:

                # calculate overlap
                overlap = calc_seq_overlap(int(reg1[START]), int(reg1[END]),
                                           int(reg2[START]), int(reg2[END]))

                # check for a an overlap
                if overlap >= OVERLAP:

                    # at this point check the evalues and build the list for
                    # the non significant regions

                    if float(reg1[EVAL]) <= float(reg2[EVAL]):

                        if ((reg2[RFAM_ACC], reg2[SEQ_ACC],
                             reg2[START]) not in non_sig_regs):
                            non_sig_regs.append(
                                (reg2[RFAM_ACC], reg2[SEQ_ACC], reg2[START]))
                    else:

                        if ((reg1[RFAM_ACC], reg1[SEQ_ACC],
                             reg1[START]) not in non_sig_regs):
                            non_sig_regs.append(
                                (reg1[RFAM_ACC], reg1[SEQ_ACC], reg1[START]))

                elif overlap is None:
                    log.debug("reg1: %s" % '\t'.join(reg1))
                    log.debug("reg2: %s" % '\t'.join(reg2))

        index = index + 1
        comp_regs = None

    return non_sig_regs

# -----------------------------------------------------------------------------


def complete_clan_seqs(sorted_clan):
    """
    Parses a sorted clan file and generates a list of regions per rfam_acc,
    which are then competed by compete_seq_regions

    sorted_clan: A valid path to a sorted clan file
    """

    fp = open(sorted_clan, 'r')

    # log regions in which calculate overlap returns None
    logging.basicConfig(
        filename="missed_overlaps.log", filemode='w', level=logging.DEBUG)

    non_sig_regs = []
    regions = []

    # read first 2 regions
    seq_prev = fp.readline().strip().split('\t')
    seq_next = fp.readline().strip().split('\t')

    # read while there are no duplicates
    while len(seq_next) > 1:

        # the same accession - create region list, otherwise the sequence is
        # significant...
        while(len(seq_next) > 1 and string.find(seq_prev[SEQ_ACC],
                                                seq_next[SEQ_ACC]) != -1):
            # add the previous only in the first occurence

            if len(regions) == 0:
                regions.append(seq_prev)  # add the first
                regions.append(seq_next)
            # else only add the next one
            else:
                regions.append(seq_next)  # add the next

            seq_prev = seq_next
            seq_next = fp.readline().strip().split('\t')

        # at some point create the seq region list
        # call compete region function
        non_sig_regs.extend(compete_seq_regions(regions, logging))

        seq_prev = seq_next
        seq_next = fp.readline().strip().split('\t')

        regions = []

    fp.close()

    # at this point update full_region table
    if len(non_sig_regs) != 0:
        db_utils.set_is_singificant_to_zero_multi(non_sig_regs)

    return non_sig_regs

# -----------------------------------------------------------------------------


def usage():
    """
    Displays information on how to run clan competition
    """

    print "\nUsage:\n------"

    print "\nclan_competition.py [clan_file|clan_dir] [-r]"

    print "\nclan_dir: A directory of sorted clan region files"
    print "clan_file: The path to a sorted clan region file"
    print "\n-r option to reset is_significant field"


# -----------------------------------------------------------------------------

if __name__ == '__main__':

    # Input will be a directory of sorted clan files or a single sorted clan
    # file

    clan_source = sys.argv[1]

    # minor input checks
    if not os.path.isdir(clan_source) and not os.path.isfile(clan_source):
        usage()
        sys.exit()

    # with -r option reset all is_significant fields back to 1
    if sys.argv.count("-r") == 1:
        print "\nReseting is_significant fields ..."
        db_utils.reset_is_significant()

    t_start = timeit.default_timer()

    print "\nCompeting Clans ...\n"

    if os.path.isdir(clan_source):

        clan_files = filter(
            lambda x: x.find(".txt") != -1, os.listdir(clan_source))

        non_sig_seqs = None

        for clan in clan_files:
            c_file = os.path.join(clan_source, clan)
            non_sig_seqs = complete_clan_seqs(c_file)

            print "%s : %s" % (str(clan[0:8]), len(non_sig_seqs))
            non_sig_seqs = None

        elapsed_time = timeit.default_timer() - t_start
        print "elapsed time: ", elapsed_time

    elif os.path.isfile(clan_source):

        c_file = clan_source
        non_sig_seqs = complete_clan_seqs(c_file)

        print "%s : %s" % (os.path.basename(clan_source).partition(".")[0],
                           len(non_sig_seqs))
        non_sig_seqs = None

        elapsed_time = timeit.default_timer() - t_start
        print "elapsed time: ", elapsed_time

    else:
        usage()
