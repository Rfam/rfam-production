#!/usr/bin/env python3
import os
import sys
import requests
import subprocess
import argparse
import hashlib

import xml.etree.ElementTree as ET

DB_RNA_TYPES = {"mirbase": "precursor RNA"}
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


def map_rnacentral_urs_wirh_db_accessions(db_accession, expert_db):
    """
    Maps a database accession with a URS accession assigned by
    RNAcentral. The limitation

    db_accession: A valid member database accession already imported
    to RNAcentral
    expert_db: RNAcentral expert database to map the SEED accessions to

    return: The corresponding RNAcentral accession (URS)
    """

    rnacentral_url = "http://wwwdev.ebi.ac.uk/ebisearch/ws/rest/rnacentral?query=\"%s\" AND expert_db:\"%s\" AND rna_type:\"%s\""

    response = requests.get(rnacentral_url % (db_accession, expert_db, DB_RNA_TYPES[expert_db.lower()]))

    rnacentral_id = None

    if response.status_code == 200:
        xml_root = ET.fromstring(response.text)
        hit_count = int(xml_root.find("hitCount").text)

        if hit_count == 1:
            rnacentral_id = xml_root.find("entries").find("entry").get("id")

    return rnacentral_id

# ---------------------------------------------------------------


def fetch_sequence_from_rnacentral(rnacentral_id):
    """
    Uses RNAcentral's API to fetch corresponding sequence based on
    RNAcentral URS id

    rnacentral_id: A valid RNAcentral URS identifier e.g.

    return: Corresponding sequence if available, none otherwise
    """

    # isolate URS if necessary
    if rnacentral_id.find('_') != -1:
        rnacentral_id = rnacentral_id.partition('_')[0]

    rnacentral_url = "https://rnacentral.org/api/v1/rna/%s.fasta"

    response = requests.get(rnacentral_url % rnacentral_id)

    sequence = None

    # if request status is OK - remove header and merge sequence
    # segments
    if response.status_code == 200:
        sequence = ''.join(response.text.strip().split('\n')[1:])

    return sequence

# ---------------------------------------------------------------

def relabel_seeds_from_rnacentral_md5_mapping(seed, dest_dir=None):
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
            miRNA_seq = miRNA_seq.replace('-', '').replace('T', 'U').replace('t', 'u')

            sequence_label = generate_seed_id_from_RNAcentral(miRNA_seq.upper())

            new_line = "\t".join([sequence_label, miRNA_seq.upper(), '\n'])

        else:
            new_line = line

        new_seed_fp.write(new_line)

    seed_fp.close()

    return new_seed_loc

# ---------------------------------------------------------------


def relabel_seeds_from_rnacentral_urs_mapping(seed, expert_db=None, dest_dir=None, clean=False):
    """
    Relabels the accessions of a SEED alignment using RNAcentral
    identifiers. This is done by matching the seed sequences, with
    sequences existing in RNAcentral using md5 hashing.

    seed: A reformatted seed in Pfam format
    expert_db: An existing RNAcentral expert database
    dest_dir: The path to the destination directory. None by default
    clean: Specifies whether to clean seeds by ommitting sequences not found in RNAcentral

    return: The path to the relabelled SEED alignement
    """

    if dest_dir is None:
        # fetch path of seed alignment
        dest_dir = os.path.split(seed)[0]

    sequence_label = 0
    new_line = ''

    write_log = False

    filename = os.path.split(seed)[1].partition('.')[0]

    new_seed_loc = os.path.join(dest_dir, filename + '_relabelled')
    seed_fp = open(seed, 'r')
    new_seed_fp = open(new_seed_loc, 'w')
    log_fp = None

    seed_seq_id = ''

    unique_seed_accs = {}

    for line in seed_fp:
        # check if this is an actual sequence line
        if line[0] != '#' and len(line) > 1 and line[0:2] != '//':
            line_elements = [x for x in line.strip().split(' ') if x != '']

            # replace alignment characters
            seed_seq_id = line_elements[0].split('/')[0]
            seed_seq = line_elements[1].replace('.', '')
            seed_seq = seed_seq.replace('-', '').replace('T', 'U').replace('t', 'u').upper()

            rnacentral_id = map_rnacentral_urs_wirh_db_accessions(seed_seq_id, expert_db)

            # if any of the sequences isn't found prepare to log missing accessions
            # sets write_log to only open the file once
            if rnacentral_id is None and write_log is False:
                seed_filename = os.path.basename(seed).partition('.')[0]
                log_fp = open(os.path.join(dest_dir, seed_filename)+'.log', 'w')
                write_log = True

            # write log file if one exists and jump to next iteration
            if rnacentral_id is None:
                log_fp.write("RNACENTRAL ID MISMATCH: %s\n" % seed_seq_id)
                continue

            rnacentral_sequence = fetch_sequence_from_rnacentral(rnacentral_id)
            coordinates = fetch_seed_sequence_coordinates(seed_seq, rnacentral_sequence)

            new_label = ''

            # make sure subsequence was found
            if (coordinates[0] != 0 or coordinates[1] != 0):
                new_label = rnacentral_id + '/' + str(coordinates[0]) + '-' + str(coordinates[1])

                if new_label not in unique_seed_accs:
                    unique_seed_accs[new_label] = ''

                # if we found a duplicated id, write log and skip this sequence
                else:
                    if write_log is False:
                        seed_filename = os.path.basename(seed).partition('.')[0]
                        log_fp = open(os.path.join(dest_dir, seed_filename) + '.log', 'w')
                        write_log = True

                    log_fp.write("RNACENTRAL DUPLICATED ID: %s\n" % new_label)
                    continue

            # if we reached this point, it means there was an id match, but no exact sequence match
            else:

                # try mapping SEED sequences in smaller segments
                if write_log is False:
                    seed_filename = os.path.basename(seed).partition('.')[0]
                    log_fp = open(os.path.join(dest_dir, seed_filename) + '.log', 'w')
                    write_log = True

                log_fp.write("RNACENTRAL SEQ MISMATCH: %s\n" % seed_seq_id)

                continue
                #else:
                    # here we need to generate the new label with coords and
                #   print (sequence)

            new_line = "\t".join([new_label, line_elements[1], '\n'])

        else:
            new_line = line

        new_seed_fp.write(new_line)

    seed_fp.close()
    new_seed_fp.close()

    # close log file if one exists
    if write_log is True:
        log_fp.close()

    return new_seed_loc

# ---------------------------------------------------------------


def map_sequence_segments(seed_seq, rnac_seq, no_segments=4):
    """
    Splits a seed sequence into smaller segments and maps the
    individual segments to the RNAcentral sequence
    param seed_seq:
    param rnac_seq:
    no_segments:

    return: RNAcentral subsequence if segments match by 75%, None otherwise
    """

    segment_hits = {}
    output = {}
    seq_match_score = 0
    seed_length = len(seed_seq)

    reference_start = None

    # calculate number of remaining nts the last segment will be extended by
    remainder = seed_length % no_segments

    # calculate even segment sizes based on no_segments
    segment_size = int((seed_length - remainder) / no_segments)


    # initialize all segment hits to 0
    index = 0
    start = 0
    end = segment_size

    # initialize
    while index < no_segments:
        # check if segment is not the last
        if index != no_segments:
            segment_hits[index + 1] = [seed_seq[start:end], 0]
        else:
            segment_hits[index + 1] = [seed_seq[start:], 0]

        start = end
        end = end + segment_size

        index += 1

    # Now map the segments to the original sequence
    for position in sorted(segment_hits.keys()):
        segment = segment_hits[position][0]

        coords = fetch_seed_sequence_coordinates(segment, rnac_seq)

        # set mapping score to 1 - the segment matches
        if coords[0] != 0 or coords[1] != 0:
            segment_hits[position][1] = 1
            seq_match_score += 1

        # set reference_start
        #if position == 1:
            #reference_start = coords[0]

    # get percentage of matches
    percentage = seq_match_score * no_segments / 100

    # check if at least 75% of the sequence segments
    # match the target sequence
    if percentage >= 75:
        """
        # need to split this to more cases
        return rnac_seq[reference_start: reference_start + seed_length]
        """
        pass

    return None

# ---------------------------------------------------------------


def build_temporary_cm_from_seed(seed_file, dest_dir=None):
    """
    Build a temporary covariance model based on the seed_file,
    which is provided as input using cmbuild.

    seed_file: Seed alignment in Stockhold format to be used as
    input to cmbuild

    dest_dir: Destination directory where the output will be
    generated

    return: True if the covariance model exists, False otherwise
    """

    if dest_dir is None:
        dest_dir = os.path.split(seed_file)[0]

    filename = os.path.basename(seed_file).partition('.')[0]

    cm_file = os.path.join(dest_dir, filename+'.cm')

    cmd = "cmbuild -F %s %s" % (cm_file, seed_file)

    subprocess.call(cmd, shell=True)

    if os.path.exists(cm_file):
        return True

    return False

# ---------------------------------------------------------------

def parse_arguments():
    """
    Basic argument parsing

    return: Argparse parser object
    """

    parser = argparse.ArgumentParser(description='Script to relabel SEED alignments')

    required_arguments = parser.add_argument_group("required arguments")
    required_arguments.add_argument("--seed", help="SEED alignment in stockholm format to relabel", type=str)

    # mutually exclusive arguments
    mutually_exclusive_arguments = parser.add_mutually_exclusive_group()
    mutually_exclusive_arguments.add_argument("--seqdb",
                                              help="Sequence file in fasta format from where to extract coordinates",
                                              type=str)
    #mutually_exclusive_arguments.add_argument("--rnac",
    #                                         help="Sets RNAcentral as the source sequence database", action="store_true")

    # group together related arguments
    rnac_option_group = parser.add_argument_group("rnacentral options")
    rnac_option_group.add_argument("--rnac", help="Sets RNAcentral as the source sequence database",
                                   action="store_true")
    rnac_option_group.add_argument("--expert-db", help="Name of experct RNAcentral database",
                                   action="store")
    rnac_option_group.add_argument("--md5", help="Map sequences based on md5 hashes",
                                   action="store_true")
    rnac_option_group.add_argument("--clean-seed", help="Ommit any SEED sequences not matching RNAcentral ones",
                                   action="store_true")

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
        if args.md5:
            reformatted_pfam_seed = relabel_seeds_from_rnacentral_md5_mapping(new_pfam_seed)
        else:
            reformatted_pfam_seed = relabel_seeds_from_rnacentral_urs_mapping(new_pfam_seed, args.expert_db,
                                                                              dest_dir=dest_dir, clean=args.clean_seed)

    # reformat to stockholm
    reformatted_stk = pfam_to_stockholm_format(reformatted_pfam_seed, dest_dir=dest_dir)

    if reformatted_stk is None:
        sys.exit("\nReformatted stockholm could not be generated!")