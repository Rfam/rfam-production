#!/usr/bin/python
'''
Created on 24 Feb 2016

@author: ikalvari
'''

import os
import json
import subprocess
import urllib2
import logging
import sys
import re
import string

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

ESL_LOCAL = "/Users/ikalvari/Downloads/infernal-1.1/easel/miniapps/esl-sfetch"
ESL_FSEQ_PATH = "/nfs/production/xfam/rfam/software/bin/esl-sfetch"
TMP_PATH = "/tmp"
FSR_PATH = "/nfs/research2/nobackup/rfam/exec/faSomeRecords"
FSR_LOCAL = "/Users/ikalvari/Desktop/faSomeRecords"
ENA_URL = "http://www.ebi.ac.uk/ena/data/view/%s&display=fasta&range=%s-%s"

ESL_PATH = None

if LSF_MODE is False:
    ESL_PATH = ESL_LOCAL
elif LSF_MODE is True:
    ESL_PATH = ESL_FSEQ_PATH
else:
    sys.exit("\nLSF_MODE has not been set properly.")

# -----------------------------------------------------------------------------


def rnac_to_json(rfam2rnac_file, fasta_dir, no_seqs=None, out_dir=None):
    '''
        This was initially developed for processing the entire Rfam2RNAcentral
        export with the output split to multiple output files with the number
        of sequences per file set by the parameter no_seqs.

        rfam2rnac_file:  Rfam2RNAcentral db dump
        fasta_dir:       The path to the directory containing the fasta files
                         of the current Rfam release
        no_seqs:         The number of sequences to split input file to
        out_dir:         The path to the output directory

    '''

    json_obj_list = []
    sequence = None

    # open a log file for tracking the obsolete sequences
    logging.basicConfig(
        filename='empty_seqs.log', filemode='w', level=logging.DEBUG)

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

            cmd = '%s %s %s' % (ESL_PATH, fam_fa_path, seq_id)
            proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)

            seq = proc.communicate()[0]

            # get sequence
            sequence = ''
            seq_bits = seq.split('\n')[1:]
            sequence = sequence.join(seq_bits)

        # if no fasta file found, get the sequence from ENA
        if (sequence == ''):
            sequence = fetch_seq_from_ena(entry)

        # build dictionary and update sequence list
        if (sequence != '' and seq_validator(sequence) is True):
            json_obj_list.append(build_json_dict(entry, sequence))
        else:
            # log obsolete sequences
            logging.debug("%s", "\t".join(entry))

        if (len(json_obj_list) == no_seqs):
            fp_out = open(
                os.path.join(out_dir, filename + str(f_index) + ".json"), 'w')

            json.dump(
                json_obj_list, fp_out, indent=2, separators=(',', ':'))

            fp_out.close()

            f_index = f_index + 1
            json_obj_list = []

        cmd = ""

    # if list != empty write remaining seqs to new file
    if (len(json_obj_list) > 0):
        fp_out = open(
            os.path.join(out_dir, filename + str(f_index) + ".json"), 'w')
        json.dump(
            json_obj_list, fp_out, indent=2, separators=(',', ':'))

        fp_out.close()

    rnac_fp.close()

# -----------------------------------------------------------------------------


def rnac_to_json_multi(seq_dir, fasta_dir, out_dir=None):
    '''
        This is an implementation of the rnac_to_json function with the
        difference that input is split to smaller files prior to the json
        generation. It exports the sequences out of the Rfam's currenct version
        of fasta files.

        seq_dir:    The path to the directory containing multiple sequence
                    files to be converted to json
        fasta_dir:  The path to the directory containing the fasta files of the
                    current Rfam release
        out_dir:    The path to the output directory
    '''

    if(out_dir is None):
        out_dir = seq_dir

    seq_files = os.listdir(seq_dir)
    seq_files = filter(lambda x: string.find(x, ".out") != -1, seq_files)

    # open a new log file and keep track of the sequences not found in the
    logging.basicConfig(filename=os.path.join(out_dir, 'obsolete_seqs.log'),
                        filemode='w', level=logging.DEBUG)

    json_obj_list = []

    for seq_file in seq_files:

        fp = open(os.path.join(seq_dir, seq_file), 'r')

        for entry in fp:
            entry = entry.strip().split('\t')

            fam_fa_path = os.path.join(fasta_dir, entry[RFAM_ACC] + ".fa")
            sequence = ''
            seq_bits = None

            # check if the fasta file exists
            if os.path.exists(fam_fa_path):
                seq_id = entry[SEQACC] + '/' + \
                    entry[SEQ_START] + '-' + entry[SEQ_END]

                cmd = '%s %s %s' % (ESL_PATH, fam_fa_path, seq_id)

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
            if (sequence != '' and seq_validator(sequence) is True):
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
    '''
        RNAcentral specific method to build the json dictionary for each entry.
        Sequences are provided as a parameter as they are exported using
        esl-sfetch and ENA via the url API.

        entry:    A list of the fields in a DB entry resulting from
                  Rfam2RNAcentral export
        sequence: Entry's corresponding sequence
    '''

    edict = {}
    species = ""

    edict["parent_accession"] = entry[SEQACC].partition('.')[0]
    edict["seq_version"] = entry[VERSION]
    edict["feature_location_start"] = entry[SEQ_START]
    edict["feature_location_end"] = entry[SEQ_END]

    edict["ncrna_class"] = entry[NCRNA_CLASS]
    edict["feature_type"] = entry[RNA_TYPE]

    if (entry[ALIGNMENT] == "seed"):
        edict["is_seed"] = '1'
    else:
        edict["is_seed"] = '0'

    edict["primary_id"] = entry[RFAM_ACC]
    edict["optional_id"] = entry[RFAM_ID]
    edict["description"] = entry[DESCRIPTION]

    edict["sequence"] = sequence

    edict["lineage"] = entry[TAX_STR].replace(';', '')

    species = entry[SPECIES]
    common_name = ""

    if(species.find('(') != -1):
        if (species.count('(') > 1):
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
    edict["ontology"] = ontology

    return edict

# -----------------------------------------------------------------------------


def fetch_seq_from_ena(entry):
    '''
        This function uses the URL API of the European Nucleotide Archive to
        fetch sequence regions according to the provided rfamseq_acc in entry.

        entry: A list of the fields in a DB entry resulting from
               Rfam2RNAcentral.pl export
    '''

    # check the strand
    start = int(entry[SEQ_START])
    end = int(entry[SEQ_END])

    seq_acc = entry[SEQACC].partition('.')[0]

    if (start < end):
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
    '''
        Checks if the sequence provided is valid fasta sequence. Returns True
        if the sequence is valid, otherwise returns False.

        sequence: A string for validation
    '''

    # checks for ascii characters that should not appear in a fasta sequence
    seq_val = re.compile(r"[.-@|\s| -)|z-~|Z-`|EFIJLOPQX|efijlopqx+,]+")

    if(seq_val.search(sequence) is None):
        return True

    False

# -----------------------------------------------------------------------------


def fa_some_records_to_json(seq_dir, fasta_dir, out_dir=None):
    '''
        This is a slightly different version of the rnac_to_json methods,
        calling UCSCs faSomeRecords executable to retrieve sequences out of
        fasta input files.

        seq_dir:   The path to the directory containing multiple sequence files
                   to be converted to json
        fasta_dir: The path to the directory containing the fasta files of the
                   current Rfam release
        out_dir:   The path to the output directory
    '''

    if(out_dir is None):
        out_dir = seq_dir

    seq_files = os.listdir(seq_dir)
    seq_files = filter(lambda x: string.find(x, ".out") != -1, seq_files)

    # open a new log file for obsolete sequences
    logging.basicConfig(filename=os.path.join(out_dir, 'obsolete_seqs.log'),
                        filemode='w', level=logging.DEBUG)
    # new family seqs
    logging.basicConfig(filename=os.path.join(out_dir, 'newfamily_seqs.log'),
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

                lfile_path = os.path.join(TMP_PATH, entry[SEQACC] + '.list')
                ofile_path = os.path.join(TMP_PATH, entry[SEQACC] + '.out')
                f_temp = open(lfile_path, 'w')
                f_temp.write(seq_id)
                f_temp.close()

                cmd = '%s %s %s %s' % (
                    FSR_LOCAL, fam_fa_path, lfile_path, ofile_path)

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
                if(sequence != '' and seq_validator(sequence) is True):
                    json_obj_list.append(build_json_dict(entry, sequence))
                else:
                    # log obsolete sequence
                    logging.debug("%s", "\t".join(entry))

            else:
                # update new families log
                logging.warning("%s", "\t".join(entry))

        fp_out = open(
            os.path.join(out_dir, seq_file.partition(".")[0] + ".json"), 'w')

        json.dump(json_obj_list, fp_out, indent=2, separators=(',', ':'))

        json_obj_list = []

        fp_out.close()
        fp.close()


# -----------------------------------------------------------------------------

def usage():
    '''
        Prints out guidelines on how to run rnac2json script.
    '''

    print "\nUsage:\n-----"

    print "\nrnac2json.py seq_dir fasta_files [out_dir]"

    print "\nseq_dir: A directory of multiple (Rfam2RNAcentral.pl) db \
            dump files"

    print "fasta_files: The directory of fasta_files of the current RFAM release"

    print "\nrnac2json.py -h for help\n"


# -----------------------------------------------------------------------------

if __name__ == '__main__':

    if (sys.argv[1] == '-h'):
        usage()

    elif(len(sys.argv) >= 3):
        seq_dir = sys.argv[1]    # directory of rfam2rnac dump files
        fasta_dir = sys.argv[2]  # directory of fasta_files
        out_dir = None           # output directory

        if (len(sys.argv) == 4):
            out_dir = sys.argv[3]

        rnac_to_json_multi(seq_dir, fasta_dir, out_dir)

    else:
        usage()

    # NOTE: when the out_dir is None the output is generated within the input
    # directory
