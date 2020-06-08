import os
import sys
import requests
import subprocess
import argparse
import hashlib

# ---------------------------------------------------------------


def fetch_seed_sequence_coordinates(seed_seq, full_seq):
    """
    Returns seed sequence start and end coordinates based on the
    accession provided as input

    accession: A valid GenBank|ENA\RNAcentral accession

    return: A tuple with start and end coordinates in this order
    """
    start = 0
    end = 0

    # use accession to fetch the coordinates - sequence might be needed
    start = full_seq.find(seed_seq)

    end = start + len(seed_seq)

    # return start and end coordinates if subsequence was found
    if start != -1:
        return (start, end)

    return (0, 0)

# ---------------------------------------------------------------


def load_fasta_file_to_dict(fasta):
    """
    Loads a fasta file (seqdb) into a dictionary with the sequence
    accession used as a key and the sequence as a value

    fasta: A valid sequence file in fasta format

    return: A python dictionary with accession:sequence pairs
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
            # replace Us in sequences (SEED) with Ts
            line = line.strip().replace('U', 'T')
            fasta_dict[accession] += line

    fasta_fp.close()

    return fasta_dict

# ---------------------------------------------------------------


def stockhom_to_pfam_format(stk_msa, dest_dir=None):
    """
    Converts a stockholm MSA to the Pfam format

    stk_msa: A valid MSA in strockholm format
    dest_dir: The destination directory where the new MSA
    will be generated

    return: The output MSA in Pfam format, None otherwise
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

    pfam_msa: A valid MSA in Pfam format
    dest_dir: The destination directory where the new MSA
    will be generated

    return: The output MSA in Stockholm format, None otherwise
    """

    filename = os.path.basename(pfam_msa).partition('.')[0]
    output_stk_msa = os.path.join(dest_dir, filename + ".stk")

    cmd = "esl-reformat stockholm %s > %s" % (pfam_msa, output_stk_msa)

    subprocess.call(cmd, shell=True)

    if not os.path.exists(output_stk_msa):
        return None

    return output_stk_msa

# ---------------------------------------------------------------


def relabel_seed_accessions(seed, accession_coords, dest_dir = None):
    """
    Re-writes a seed file with the sequence coordinates appended to the
    sequence accessions

    seed: A seed file with no star-end sequence coordinates

    return: The path to the newly relabelled SEED alignment
    """

    if dest_dir is None:
        # fetch path of seed alignment
        dest_dir = os.path.split(seed)[0]

    sequence_label = 0
    new_line = ''

    filename = os.path.split(seed)[1].partition('.')[0]

    new_seed_loc = os.path.join(dest_dir, filename+'_relabelled')
    seed_fp = open(seed, 'r')
    new_seed_fp = open(new_seed_loc, 'w')

    for line in seed_fp:
        # check if this is an actual sequence line
        if line[0] != '#' and len(line) > 1 and line[0:2] != '//':
            line_elements = [x for x in line.strip().split(' ') if x != '']

            sequence_label = line_elements[0] + '/' + accession_coords[line_elements[0]]

            new_line = "\t".join([sequence_label, line_elements[1], '\n'])
        else:
            new_line = line

        new_seed_fp.write(new_line)

    seed_fp.close()

    return new_seed_loc

# ---------------------------------------------------------------


def fetch_RNAcentral_id(sequence):
    """
    Looks for a sequence match in RNAcentral based on sequence md5
    and fetches the corresponding RNAcentral accession

    sequence: A valid DNA/RNA sequence

    return: Returns RNAcentral id, otherwise returns None
    """

    sequence_md5 = sequence_to_md5(sequence)
    rnacentral_url = 'https://rnacentral.org/api/v1/rna'
    response = requests.get(rnacentral_url, params={'md5': sequence_md5})

    data = response.json()

    if data['count'] > 0:
        return data['results'][0]['rnacentral_id']

    return None

# ---------------------------------------------------------------


def generate_seed_id_from_RNAcentral(sequence):
    """
    Generates a seed accession based on a sequence mad5 match in RNAcentral

    sequence: A valid DNA/RNA sequence

    return: Returns RNAcentral id, otherwise returns None
    """

    sequence_md5 = sequence_to_md5(sequence)
    rnacentral_url = 'https://rnacentral.org/api/v1/rna'
    response = requests.get(rnacentral_url, params={'md5': sequence_md5})

    data = response.json()

    if data['count'] > 0:
        return data['results'][0]['rnacentral_id'] + "/1-" + str(data['results'][0]['length'])

    return None

# ---------------------------------------------------------------


def sequence_to_md5(sequence):
    """
    Converts a sequence to an md5 hash after replacing Us with
    Ts

    sequence: A valid RNA/DNA sequence

    return: MD5 hash of the sequence
    """

    md5_converter = hashlib.md5()
    # convert to DNA
    sequence = sequence.replace('U', 'T')
    md5_converter.update(sequence.encode('utf-8'))
    sequence_md5 = md5_converter.hexdigest()

    return sequence_md5

# ---------------------------------------------------------------


def validate_sequences(seed_sequence, extracted_full_seq):
    """
    Validates whether the SEED sequence matches the sequence
    extracted at specific coordinates

    seed_sequence: A DNA/RNA sequecne extracted from the SEED alignment
    extracted_full: A DNA/RNA subsequence extracted at specific locations

    return: True if the sequences match, False otherwise. Returns False by default
    """

    new_seed_sequence = seed_sequence.replace('U', 'T')
    if extracted_full_seq.find(new_seed_sequence) != -1:
        return True

    return False

# ---------------------------------------------------------------


def seed_to_fasta(seed_msa, dest_dir=None):
    """
    Converts a multiple sequence alignment (MSA) to fasta

    param seed_msa: A valid Rfam SEED file in stockholm format to convert to fasta

    return: Path to updated seed file
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


def relabel_seeds_from_rnacentral(seed, dest_dir=None):
    """
    Relabels the accessions of a SEED alignment using RNAcentral
    identifiers. This is done by matching the seed sequences, with
    sequences existing in RNAcentral using md5 hashing.

    seed: A reformatted seed in Pfam format
    dest_dir: The path to the destination directory. None by default

    return: The path to the relabelled SEED alignement
    """

    if dest_dir is None:
        # fetch path of seed alignment
        dest_dir = os.path.split(seed)[0]

    sequence_label = 0
    new_line = ''

    filename = os.path.split(seed)[1].partition('.')[0]

    new_seed_loc = os.path.join(dest_dir, filename+'_relabelled')
    seed_fp = open(seed, 'r')
    new_seed_fp = open(new_seed_loc, 'w')

    for line in seed_fp:
        # check if this is an actual sequence line
        if line[0] != '#' and len(line) > 1 and line[0:2] != '//':
            line_elements = [x for x in line.strip().split(' ') if x != '']

            # replace alignment characters
            miRNA_seq = line_elements[1].replace('.', '')
            miRNA_seq = line_elements[1].replace('-', '')

            sequence_label = generate_seed_id_from_RNAcentral(miRNA_seq)

            new_line = "\t".join([sequence_label, miRNA_seq, '\n'])
        else:
            new_line = line

        new_seed_fp.write(new_line)

    seed_fp.close()

    return new_seed_loc

# ---------------------------------------------------------------


def parse_arguments():
    """
    Basic argument parsing

    return: Argparse parser object
    """

    parser = argparse.ArgumentParser(description='Script to relabel SEED alignments')

    parser.add_argument("--seed", help="SEED alignment in stockholm format to relabel", type=str)

    mutually_exclusive_arguments = parser.add_mutually_exclusive_group()
    mutually_exclusive_arguments.add_argument("--seqdb",
                                              help="Sequence file in fasta format from where to extract coordinates", type=str)
    mutually_exclusive_arguments.add_argument("--rnac",
                                              help="Sets RNAcentral as the source sequence database", action="store_true")

    return parser

# ---------------------------------------------------------------

if __name__ == '__main__':

    parser = parse_arguments()
    args = parser.parse_args()

    # temporarily use seed source dir
    dest_dir = os.path.split(args.seed)[0]

    # declaring a variable to store the path to the reformatted SEED alignment
    reformatted_pfam_seed = None

    # convert seed to fasta
    seed_fasta = seed_to_fasta(args.seed, dest_dir=None)

    # convert stockholm to pfam
    new_pfam_seed = stockhom_to_pfam_format(args.seed, dest_dir=dest_dir)

    if not args.rnac:
        # load sequences to dict
        seed_seq_dict = load_fasta_file_to_dict(seed_fasta)
        full_seq_dict = load_fasta_file_to_dict(args.seqdb)

        # construct accession coords dictionary
        accession_coords = {}
        for accession in seed_seq_dict.keys():

            (start, end) = fetch_seed_sequence_coordinates(seed_seq_dict[accession], full_seq_dict[accession])

            # validate sequences
            check = validate_sequences(seed_seq_dict[accession], full_seq_dict[accession][start:end])

            # check if coordinates were extracted successfully,
            if end != 0 and check is True:
                # start point is one position shifted to the right on actual sequence
                start += 1
                accession_coords[accession] = '-'.join((str(start), str(end)))
            else:
                print ("Unable to extract coordinates for accession: %s" % accession)


        # now rewrite the seed labels
        reformatted_pfam_seed = relabel_seed_accessions(new_pfam_seed, accession_coords, dest_dir=dest_dir)

    # relabel SEED accessions using RNAcentral identifiers
    else:
        reformatted_pfam_seed = relabel_seeds_from_rnacentral(new_pfam_seed)


    # reformat to stockholm
    reformatted_stk = pfam_to_stockholm_format(reformatted_pfam_seed, dest_dir=dest_dir)

    if reformatted_stk is None:
        sys.exit("\nReformatted stockholm could not be generated!")
