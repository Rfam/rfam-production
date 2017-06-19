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
Description:    Script to generate fasta files for all family regions
                in full_region. For execution on lsf call fasta_gen_handler
                that generates distinct shell scripts per family to enable
                recovery.
"""

# ---------------------------------IMPORTS-------------------------------------

import os
import sys
import subprocess
import logging
import re
import gzip
from utils import RfamDB
from config import rfam_config

# -----------------------------------------------------------------------------
# To be modified accordingly
LSF = False

RFAM_ACC = 0
SEQ_ACC = 1
START = 2
END = 3
DESC = 4

ESL_PATH = None

ESL_LSF = rfam_config.ESL_PATH
ESL_LOCAL = rfam_config.ESL_PATH

if LSF is True:
    ESL_PATH = ESL_LSF
else:
    ESL_PATH = ESL_LOCAL


# -----------------------------------------------------------------------------

def generate_fasta(seq_file, out_dir):
    """
    Uses esl-sfetch to generate family specific fasta files out of seq_file
    which is provided as source (e.g. rfamseq11.fa). It will generate fasta
    files for all families by default

    seq_file:   The path to rfamseq input file in fasta format, for
                generating the fasta files

    out_dir:    Destination directory where the files will be
                generated
    """

    sequence = ''
    fp_out = None
    seq_bits = None

    # logging sequences not exported
    # rename this to family log
    log_file = os.path.join(out_dir, "missing_seqs.log")
    logging.basicConfig(
        filename=log_file, filemode='w', level=logging.INFO)

    # connect to db
    cnx = RfamDB.connect()

    # get a new buffered cursor
    cursor = cnx.cursor(raw=True)

    # fetch clan specific family full_region data and sequence description

    query = ("SELECT fr.rfam_acc, fr.rfamseq_acc, fr.seq_start, fr.seq_end, rf.description\n"
             "FROM full_region fr, rfamseq rf\n"
             "WHERE fr.rfamseq_acc=rf.rfamseq_acc\n"
             "AND fr.is_significant=1\n"
             "ORDER BY fr.rfam_acc")

    # execute the query
    cursor.execute(query)

    for region in cursor:

        # new family
        if str(region[RFAM_ACC]) != rfam_acc:
            # check if there's no open file
            if fp_out is not None:
                fp_out.close()

            # open new fasta file
            fp_out = gzip.open(
                os.path.join(out_dir, str(region[RFAM_ACC]) + ".fa.gz"), 'w')

        rfam_acc = region[RFAM_ACC]

        cmd = "%s -c %s/%s %s %s" % (ESL_PATH,
                                     str(region[START]), str(region[END]),
                                     seq_file, str(region[SEQ_ACC]))

        proc = subprocess.Popen(
            cmd, shell=True, stdout=subprocess.PIPE)

        seq = proc.communicate()[0]

        # get sequence
        sequence = ''
        seq_bits = seq.split('\n')[1:]
        sequence = sequence.join(seq_bits)

        # print sequence

        if sequence != '' and seq_validator(sequence) is True:
            # write header
            fp_out.write(">%s/%s-%s %s\n" % (str(region[SEQ_ACC]),
                                             str(region[START]),
                                             str(region[END]),
                                             str(region[DESC])))

            # write sequence
            fp_out.write(sequence + '\n')

        else:
            # logging sequences that have not been exported
            logging.info(sequence)

    # close last file
    fp_out.close()

    # disconnect from DB
    cursor.close()
    RfamDB.disconnect(cnx)

# -----------------------------------------------------------------------------


def generate_fasta_single(seq_file, rfam_acc, out_dir):
    """
    Uses esl-sfetch to generate family specific fasta files out of seq_file
    which is provided as source. Works on single family based on rfam_acc.
    Files are generated in a compressed .fa.gz format

    seq_file:   This is the the path to rfamseq input file in fasta format,
                for generating the fasta files

    rfam_acc:   The rfam_acc of a specific family

    out_dir:    This is the destination directory where the files will be
                generated
    """

    sequence = ''
    fp_out = None
    seq_bits = None

    # logging sequences not exported
    # rename this to family log
    log_file = os.path.join(out_dir, rfam_acc + ".log")
    logging.basicConfig(
        filename=log_file, filemode='w', level=logging.INFO)

    # connect to db
    cnx = RfamDB.connect()

    # get a new buffered cursor
    cursor = cnx.cursor(raw=True)

    # fetch sequence accessions for specific family - significant only!!
    query = ("SELECT fr.rfam_acc, fr.rfamseq_acc, fr.seq_start, fr.seq_end, rf.description\n"
             "FROM full_region fr, rfamseq rf\n"
             "WHERE fr.rfamseq_acc=rf.rfamseq_acc\n"
             "AND fr.is_significant=1\n"
             "AND fr.rfam_acc=\'%s\'") % (rfam_acc)

    # execute the query
    cursor.execute(query)

    # open a new fasta output file
    fp_out = gzip.open(
        os.path.join(out_dir, str(rfam_acc) + ".fa.gz"), 'w')

    for region in cursor:

        cmd = "%s -c %s/%s %s %s" % (ESL_PATH,
                                     str(region[START]), str(region[END]),
                                     seq_file, str(region[SEQ_ACC]))

        proc = subprocess.Popen(
            cmd, shell=True, stdout=subprocess.PIPE)

        seq = proc.communicate()[0]

        # get sequence
        sequence = ''
        seq_bits = seq.split('\n')[1:]
        sequence = sequence.join(seq_bits)

        # print sequence

        if sequence != '' and seq_validator(sequence) is True:
            # write header
            fp_out.write(">%s/%s-%s %s\n" % (str(region[SEQ_ACC]),
                                             str(region[START]),
                                             str(region[END]),
                                             str(region[DESC])))

            # write sequence
            fp_out.write(sequence + '\n')

        else:
            # logging sequences that have not been exported
            logging.info(str(region[SEQ_ACC]))

    # close last file
    fp_out.close()

    # disconnect from DB
    cursor.close()
    RfamDB.disconnect(cnx)


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

    False

# -----------------------------------------------------------------------------


def usage():
    """
    Displays information on how to run fasta_generator
    """

    print "\n Usage\n-------"
    print "\npython fasta_generator.py seq_file rfam_acc out_dir"
    print "\n-h option for help\n"


# -----------------------------------------------------------------------------

if __name__ == '__main__':

    # need to add more parameters for local and lsf execution rather that global
    # variable

    # some parameter checking
    if len(sys.argv) > 3:
        sequence_file = sys.argv[1]
        rfam_acc = sys.argv[2]
        output_dir = sys.argv[3]

        generate_fasta_single(sequence_file, rfam_acc, output_dir)

    else:
        if sys.argv[1] == "-h":
            usage()
        else:
            print "\nIncorrect Input."
            usage()
