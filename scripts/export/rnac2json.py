#!/usr/bin/python

"""
Copyright [2009-2017] EMBL-European Bioinformatics Institute
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
Requires esl-fetch to be added to PATH
"""

import json
import logging
import os
import re
import string
import subprocess
import argparse
import urllib2

from config import rfam_config as rfc

LSF_MODE = False

# DB Fields
RFAM_ACC = 0
RFAM_ID = 1
RNA_TYPE = 2
NCRNA_CLASS = 3
ALIGNMENT = 4
DESCRIPTION = 5
NCBI_ID = 6
SEQACC = 7
SEQ_START = 8
SEQ_END = 9
BITSCORE = 10
DBXREFS = 11
PMIDS = 12
VERSION = 13
SPECIES = 14
TAX_STR = 15

ENA_URL = rfc.ENA_URL
TMP_PATH = rfc.TMP_PATH

# -------------------------------------------------------------------------


def rnac_to_json(rfam2rnac_file, fasta_dir, no_seqs=None, out_dir=None):
    """
    This was initially developed for processing the entire Rfam2RNAcentral
    export with the output split to multiple output files with the number
    of sequences per file set by the parameter no_seqs.

    rfam2rnac_file:  Rfam2RNAcentral db dump
    fasta_dir:       The path to the directory containing the fasta files
                     of the current Rfam release
    no_seqs:         The number of sequences to split input file to
    out_dir:         The path to the output directory
    """

    json_obj_list = []
    sequence = None

    # open a log file for tracking the obsolete sequences
    logging.basicConfig(
        filename="empty_seqs.log", filemode='w', level=logging.DEBUG)

    rnac_fp = open(rfam2rnac_file, 'r')
    filename = os.path.basename(rfam2rnac_file).partition('.')[0]

    # drop header
    rnac_fp.readline()

    f_index = 1  # file index

    for entry in rnac_fp:

        # get entry fields
        entry = entry.strip().split('\t')

        # get family fa file path
        fam_fa_path = os.path.join(fasta_dir, entry[RFAM_ACC] + ".fa")

        sequence = ''
        seq_bits = None

        # check if the fasta file exists
        if os.path.exists(fam_fa_path):
            seq_id = entry[SEQACC] + '/' + \
                     entry[SEQ_START] + '-' + entry[SEQ_END]

            cmd = "esl-sfetch %s %s" % (fam_fa_path, seq_id)
            proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)

            seq = proc.communicate()[0]

            # get sequence
            sequence = ''
            seq_bits = seq.split('\n')[1:]
            sequence = sequence.join(seq_bits)

        # if no fasta file found, get the sequence from ENA
        if sequence == '':
            sequence = fetch_seq_from_ena(entry)

        # build dictionary and update sequence list
        if sequence != '' and seq_validator(sequence) is True:
            json_obj_list.append(build_json_dict(entry, sequence))
        else:
            # log obsolete sequences
            logging.debug("%s", "\t".join(entry))

        if len(json_obj_list) == no_seqs:
            fp_out = open(
                os.path.join(out_dir, filename + str(f_index) + ".json"), 'w')

            json.dump(
                json_obj_list, fp_out, indent=2, separators=(',', ':'))

            fp_out.close()

            f_index = f_index + 1
            json_obj_list = []

        cmd = ''

    # if list != empty write remaining seqs to new file
    if len(json_obj_list) > 0:
        fp_out = open(
            os.path.join(out_dir, filename + str(f_index) + ".json"), 'w')
        json.dump(
            json_obj_list, fp_out, indent=2, separators=(',', ':'))

        fp_out.close()

    rnac_fp.close()


# -----------------------------------------------------------------------------


def rnac_to_json_multi(seq_dir, fasta_dir, out_dir=None):
    """
    This is an implementation of the rnac_to_json function with the
    difference that input is split to smaller files prior to the json
    generation. It exports the sequences out of the Rfam's currenct version
    of fasta files.

    seq_dir:    The path to the directory containing multiple sequence
                files to be converted to json
    fasta_dir:  The path to the directory containing the fasta files of the
                current Rfam release
    out_dir:    The path to the output directory
    """

    if out_dir is None:
        out_dir = seq_dir

    seq_files = os.listdir(seq_dir)
    seq_files = filter(lambda x: string.find(x, ".txt") != -1, seq_files)

    # open a new log file and keep track of the sequences not found in the
    logging.basicConfig(filename=os.path.join(out_dir, "obsolete_seqs.log"),
                        filemode='w', level=logging.DEBUG)

    json_obj_list = []

    for seq_file in seq_files:

        fp = open(os.path.join(seq_dir, seq_file), 'r')
        seq_file_path = ''
        for entry in fp:
            entry = entry.strip().split('\t')

            if entry[ALIGNMENT] == 'seed' or entry[ALIGNMENT] == 'SEED':
                seq_file_path = os.path.join(fasta_dir, "Rfam.seed")
            else:
                seq_file_path = os.path.join(fasta_dir, entry[RFAM_ACC] + ".fa")

            sequence = ''
            seq_bits = None

            # check if the fasta file exists
            if os.path.exists(seq_file_path):
                seq_id = entry[SEQACC] + '/' + \
                         entry[SEQ_START] + '-' + entry[SEQ_END]

                cmd = "esl-sfetch %s %s" % (seq_file_path, seq_id)

                proc = subprocess.Popen(
                    cmd, shell=True, stdout=subprocess.PIPE)

                seq = proc.communicate()[0]

                # get sequence
                sequence = ''
                seq_bits = seq.split('\n')[1:]
                sequence = sequence.join(seq_bits)

            # if sequence not found try fetching the sequence from ENA
            # if (sequence == ''):
            #    sequence = fetch_seq_from_ena(entry)

            # check if seq string still empty,
            if sequence != '' and seq_validator(sequence) is True:
                json_obj_list.append(build_json_dict(entry, sequence))
            else:
                # log obsolete sequence
                logging.debug("%s", "\t".join(entry))

            seq_id = None

        fp_out = open(
            os.path.join(out_dir, seq_file.partition(".")[0] + ".json"), 'w')

        json.dump(
            json_obj_list, fp_out, indent=2, separators=(',', ':'))

        json_obj_list = []

        fp_out.close()
        fp.close()


# -----------------------------------------------------------------------------

def build_json_dict(entry, sequence):
    """
    RNAcentral specific method to build the json dictionary for each entry.
    Sequences are provided as a parameter as they are exported using
    esl-sfetch and ENA via the url API.

    entry:    A list of the fields in a DB entry resulting from
              Rfam2RNAcentral export
    sequence: Entry's corresponding sequence
    """

    edict = {}
    species = ''

    edict["parent_accession"] = entry[SEQACC].partition('.')[0]
    edict["seq_version"] = entry[VERSION]
    edict["feature_location_start"] = entry[SEQ_START]
    edict["feature_location_end"] = entry[SEQ_END]

    edict["ncrna_class"] = entry[NCRNA_CLASS]
    edict["feature_type"] = entry[RNA_TYPE]

    if entry[ALIGNMENT] == "seed":
        edict["is_seed"] = '1'
    else:
        edict["is_seed"] = '0'

    edict["primary_id"] = entry[RFAM_ACC]
    edict["optional_id"] = entry[RFAM_ID]
    edict["description"] = entry[DESCRIPTION]

    edict["sequence"] = sequence

    edict["lineage"] = entry[TAX_STR].replace(';', '')

    species = entry[SPECIES]
    common_name = ''

    if species.find('(') != -1:
        if species.count('(') > 1:
            species = species.rsplit('(', 1)
            common_name = species[1]
        else:
            species = species.partition('(')
            common_name = species[2]
        common_name = common_name.partition(')')
        common_name = common_name[0].capitalize()
        species = species[0].strip()

    edict["species"] = species
    edict["common_name"] = common_name

    edict["ncbi_tax_id"] = entry[NCBI_ID]

    references = map(lambda x: x[5:].strip(), entry[PMIDS].split(','))
    edict["references"] = references

    ontology = map(lambda x: x.strip(), entry[DBXREFS].split(','))
    ontology_dict = {}

    if len(ontology) > 0:
        for ontology_element in ontology:
            elements = ontology_element.partition(':')
            ontology_dict[elements[0]] = elements[2]
    
        edict["ontology"] = ontology_dict
    else:
        edict["ontology"] = ontology

    return edict


# -----------------------------------------------------------------------------


def fetch_seq_from_ena(entry):
    """
    This function uses the URL API of the European Nucleotide Archive to
    fetch sequence regions according to the provided rfamseq_acc in entry.

    entry: A list of the fields in a DB entry resulting from
           Rfam2RNAcentral.pl export
    """

    # check the strand
    start = int(entry[SEQ_START])
    end = int(entry[SEQ_END])

    seq_acc = entry[SEQACC].partition('.')[0]

    if start < end:
        seq_url = ENA_URL % (seq_acc, entry[SEQ_START], entry[SEQ_END])
    else:
        seq_url = ENA_URL % (seq_acc, entry[SEQ_END], entry[SEQ_START])

    seq = urllib2.urlopen(seq_url).read()

    # get sequence
    sequence = ''
    seq_bits = seq.split('\n')[1:]
    sequence = sequence.join(seq_bits)

    return sequence


# -----------------------------------------------------------------------------


def seq_validator(sequence):
    """
    Checks if the sequence provided is valid fasta sequence. Returns True
    if the sequence is valid, otherwise returns False.

    sequence: A string for validation
    """

    # checks for ascii characters that should not appear in a fasta sequence
    seq_val = re.compile(r"[.-@|\s| -)|z-~|Z-`|EFIJLOPQX|efijlopqx+,]+")

    if seq_val.search(sequence) is None:
        return True

    return False


# -----------------------------------------------------------------------------


def fa_some_records_to_json(seq_dir, fasta_dir, out_dir=None):
    """
    This is a slightly different version of the rnac_to_json methods,
    calling UCSCs faSomeRecords executable to retrieve sequences out of
    fasta input files.

    seq_dir:   The path to the directory containing multiple sequence files
               to be converted to json
    fasta_dir: The path to the directory containing the fasta files of the
               current Rfam release
    out_dir:   The path to the output directory
    """

    if out_dir is None:
        out_dir = seq_dir

    seq_files = os.listdir(seq_dir)
    seq_files = filter(lambda x: string.find(x, ".out") != -1, seq_files)

    # open a new log file for obsolete sequences
    logging.basicConfig(filename=os.path.join(out_dir, "obsolete_seqs.log"),
                        filemode='w', level=logging.DEBUG)
    # new family seqs
    logging.basicConfig(filename=os.path.join(out_dir, "newfamily_seqs.log"),
                        filemode='w', level=logging.WARNING)

    json_obj_list = []

    for seq_file in seq_files:

        fp = open(os.path.join(seq_dir, seq_file), 'r')

        for entry in fp:

            entry = entry.strip().split('\t')

            fam_fa_path = os.path.join(fasta_dir, entry[RFAM_ACC] + ".fa")
            sequence = ''

            # check if the fasta file exists
            if os.path.exists(fam_fa_path):
                seq_id = entry[SEQACC] + '/' + \
                         entry[SEQ_START] + '-' + entry[SEQ_END]

                lfile_path = os.path.join(TMP_PATH, entry[SEQACC] + ".list")
                ofile_path = os.path.join(TMP_PATH, entry[SEQACC] + ".out")
                f_temp = open(lfile_path, 'w')
                f_temp.write(seq_id)
                f_temp.close()

                cmd = "esl-sfetch %s %s %s" % (fam_fa_path, lfile_path, ofile_path)

                # call faSomeRecords to retrieve sequence
                subprocess.call(cmd, shell=True)

                # time.sleep(5)

                # open fsr out file and get sequence
                fsr_out = open(ofile_path, 'r')

                # drop header
                fsr_out.readline()

                # build seq string from file
                sequence = ''
                for line in fsr_out:
                    line = line.strip()
                    sequence = sequence + line

                fsr_out.close()

                sequence = sequence.strip()

                # remove tmp files
                os.remove(lfile_path)
                os.remove(ofile_path)

                # update json obj list
                if sequence != '' and seq_validator(sequence) is True:
                    json_obj_list.append(build_json_dict(entry, sequence))
                else:
                    # log obsolete sequence
                    logging.debug("%s", '\t'.join(entry))

            else:
                # update new families log
                logging.warning("%s", '\t'.join(entry))

        fp_out = open(
            os.path.join(out_dir, seq_file.partition('.')[0] + ".json"), 'w')

        json.dump(json_obj_list, fp_out, indent=2, separators=(',', ':'))

        json_obj_list = []

        fp_out.close()
        fp.close()


# -----------------------------------------------------------------------------

def parse_arguments():
    """
    Basic argument parsing using Python's argparse

    return: Argparse parser object
    """

    parser = argparse.ArgumentParser("Tool to convert rnacentral export to json")

    parser.add_argument("--input",
                        help="A directory of multiple (Rfam2RNAcentral.pl) dump files",
                        action="store")
    parser.add_argument("--rfam-fasta",
                        help="The directory of fasta_files of the current Rfam release",
                        action="store")
    parser.add_argument("--outdir", help='Output directory', action="store", default=None)

    return parser

# -----------------------------------------------------------------------------

if __name__ == '__main__':

    parser = parse_arguments()
    args = parser.parse_args()

    seq_dir = args.input
    fasta_dir = args.rfam_fasta
    outdir = args.outdir

    rnac_to_json_multi(seq_dir, fasta_dir, outdir)
