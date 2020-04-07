import os
import sys
import requests
import subprocess
import argparse

# ---------------------------------------------------------------


def fetch_seed_sequence_coordinates(seed_seq, full_seq):
    """
    Returns seed sequence start and end coordinates based on the
    accession provided as input

    :param accession: A valid GenBank|ENA\RNAcentral accession

    :return: A tuple with start and end coordinates in this order
    """
    start = 0
    end = 0

    # use accession to fetch the coordinates - sequence might be needed
    start = full_seq.find(seed_seq)
    end = start + len(seed_seq) - 1

    # return start and end coordinates if subsequence was found
    if start != -1:
        return (start, end)

    return (None, None)

# ---------------------------------------------------------------


def load_fasta_file_to_dict(fasta):
    """

    :param fasta:
    :return:
    """

    fasta_dict = {}

    fasta_fp = open(fasta, 'r')

    accession = ""

    flag_init_seq = 0

    for line in fasta_fp:

        # fetching accession from fasta header
        # with the 1st statement you ensure that
        if line[0] == '>':
            line = line.strip()
            # go on and extract header and initialize key pair in dict
            elements = line[1:].split(' ')
            accession = elements[0]
            fasta_dict[accession] = ""

        # keep building the sequence
        else:
            line = line.strip()
            fasta_dict[accession] += line

    fasta_fp.close()

    return fasta_dict

# ---------------------------------------------------------------


def stockhom_to_pfam_format(stk_msa, dest_dir=None):
    """
    Converts a stockholm MSA to the Pfam format

    :param stk_msa: A valid MSA in strockholm format
    :param dest_dir: The destination directory where the new MSA
    will be generated

    :return: The output MSA in Pfam format, None otherwise
    """
    filename = os.path.basename(stk_msa).partition('.')[0]
    output_pfam_msa = os.path.join(dest_dir, filename + ".pfam")

    cmd = "esl-reformat pfam %s > %s" % (stk_msa, output_pfam_msa)

    subprocess.call(cmd, shell=True)

    if not os.path.exists(output_pfam_msa):
        return None

    return output_pfam_msa

# ---------------------------------------------------------------


def pfam_to_stockholm_format(pfam_msa, dest_dir=None):
    """
    Converts a Pfam MSA to the stockholm format

    :param pfam_msa: A valid MSA in Pfam format
    :param dest_dir: The destination directory where the new MSA
    will be generated

    :return: The output MSA in Stockholm format, None otherwise
    """

    filename = os.path.basename(pfam_msa).partition('.')[0]
    output_stk_msa = os.path.join(dest_dir, filename + ".stk")

    cmd = "esl-reformat stockholm %s > %s" % (pfam_msa, output_stk_msa)

    subprocess.call(cmd, shell=True)

    if not os.path.exists(output_stk_msa):
        return None

    return output_stk_msa

# ---------------------------------------------------------------


def parse_and_rewrite_seed_alignment(seed, dest_dir = None):
    """

    :param seed:
    :return:
    """

    if dest_dir is None:
        # fetch path of seed alignment
        dest_dir = os.path.split(seed)[0]

    sequence_label = 0
    new_line = ''

    seed_fp = open(seed, 'r')
    new_seed_fp = open(os.path.join(dest_dir, "new_seed"), 'w')

    for line in seed_fp:
        # check if this is an actual sequence line
        if line[0] != '#' and len(line) > 1 and line[0:2] != '//':
            line_elements = [x for x in line.strip().split(' ') if x != '']

            (start, end) = fetch_seed_sequence_coordinates(line_elements[0])

            sequence_label = line_elements[0] + '/' + str(start) + '-' + str(end)

            new_line = "\t".join([sequence_label, line_elements[1], '\n'])
        else:
            new_line = line

        new_seed_fp.write(new_line)

    seed_fp.close()

# ---------------------------------------------------------------


def seed_to_fasta(seed_msa, dest_dir=None):
    """

    :param seed_msa:
    :return:
    """

    filename = ""

    path_elements = os.path.split(seed_msa)

    if dest_dir is None:
        dest_dir = path_elements[0]

    if "." in path_elements[1]:
        filename = path_elements[1].partition('.')[0]
    else:
        filename = path_elements[1]

    cmd = "esl-sfetch -o %s -f %s %s"

    seed_fp = open(seed_msa, 'r')
    tmp_acc_list = os.path.join("/tmp", filename + "_acc_list.txt")
    tmp_acc_fp = open(tmp_acc_list, 'w')

    for line in seed_fp:
        line = line.strip()

        if line != '' and line[0] != "#":
            accession = line.split(' ')[0]

            if accession.find('//') == -1:
                tmp_acc_fp.write(accession+'\n')

    seed_fp.close()
    tmp_acc_fp.close()

    seed_fasta_path = os.path.join(dest_dir, filename + '.fa')
    subprocess.call(cmd % (seed_fasta_path, seed_msa, tmp_acc_list), shell=True)

    # clean temp file
    os.remove(tmp_acc_list)

    # return None if the fasta file was not created
    if not os.path.exists:
        return None

    return seed_fasta_path

# ---------------------------------------------------------------


if __name__ == '__main__':

    # e.g. gammacov_3utr.fa
    seed = sys.argv[1]

    # e.g. gammacov_genomes.fa
    #all_fasta = sys.argv[2]

    #seed_seqs = load_fasta_file_to_dict(seed)
    #all_seqs = load_fasta_file_to_dict(all_fasta)

    seed_to_fasta(seed, dest_dir=None)

    """
    for accession in seed_seqs.keys():
        (start, end) = fetch_seed_sequence_coordinates(seed_seqs[accession], all_seqs[accession])

        # check if coordinates extracted
        if start is not None:
            new_pfam_msa = stockhom_to_pfam_format(seed)
            parse_and_rewrite_seed_alignment(new_pfam_msa, dest_dir=None)

    """
    #parse_and_rewrite_seed_alignment(new_pfam_msa)
