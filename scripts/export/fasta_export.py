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
Description:    Used for RT support purposes. It generates fasta files by
                running a variety of DB queries. Query must return values
                in the following order:
                rfam_acc, rfamseq_acc, seq_start, seq_end, description

TO DO:          1. Move seq_validator to utils and import it
                2. Validate mysql query
"""
# ---------------------------------IMPORTS-------------------------------------

import os
import sys
import gzip
import subprocess
import re
import logging
import argparse
from utils import RfamDB
from config import rfam_config

# ---------------------------------GLOBALS-------------------------------------

ESL_PATH = rfam_config.ESL_PATH
OUT_FILE_NAME = "seq_export.fa.gz"

RFAM_ACC = 0
SEQ_ACC = 1
START = 2
END = 3
DESC = 4

# -----------------------------------------------------------------------------


def export_sequences(seq_db, sql, filename=None, out_dir=None):
    """
    Exporting sequences from rfam_live and generating a fasta file
    by fetching the corresponding regions from seq_db provided as param

    seq_db:     A fasta sequence database to extract sequence regions from.
                Default seq_db is rfamseq11.fa
    sql:        The query to execute (string or valid .sql file)
    filename:   Ouput filename
    out_dir:    A path to the output directory
    """

    log_file = os.path.join(out_dir, "missing_seqs.log")
    logging.basicConfig(
        filename=log_file, filemode='w', level=logging.INFO)

    cnx = RfamDB.connect()
    cursor = cnx.cursor(raw=True)

    query = ''

    if os.path.isfile(sql):
        fp = open(sql, 'r')
        query = ' '.join(fp.readlines())

    else:
        query = sql

    cursor.execute(query)

    # open an output file
    fp_out = None

    if filename is not None:
        fp_out = gzip.open(
            os.path.join(out_dir, filename + ".fa.gz"), 'w')
    else:
        fp_out = gzip.open(os.path.join(out_dir, OUT_FILE_NAME), 'w')

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

        if (sequence != '' and seq_validator(sequence) is True):
            # write header
            fp_out.write(">%s/%s-%s %s\n" % (str(region[SEQ_ACC]),
                                             str(region[START]),
                                             str(region[END]),
                                             str(region[DESC])))

            # write sequence
            fp_out.write(sequence + '\n')

        else:
            logging.info(sequence)

    fp_out.close()

    cursor.close()
    RfamDB.disconnect(cnx)

# -----------------------------------------------------------------------------


def seq_validator(sequence):
    """
    Checks if the sequence provided is valid fasta sequence. Returns True
    if the sequence is valid, otherwise returns False

    sequence:   A string for validation
    """

    # checks for ascii characters that should not appear in a fasta sequence
    seq_val = re.compile(r"[.-@|\s| -)|z-~|Z-`|EFIJLOPQX|efijlopqx+,]+")

    if(seq_val.search(sequence) is None):
        return True

    return False

# -----------------------------------------------------------------------------


def usage():
    """
    Parses arguments and displays usage information on screen
    """

    parser = argparse.ArgumentParser(
        description="Rfam fasta export tool", epilog='')

    # group required arguments together
    req_args = parser.add_argument_group("required arguments")

    req_args.add_argument("--sql", help="query to execute (string or .sql file)",
                          type=str, required=True)

    parser.add_argument(
        "--infile", help="sequence database file (rfamseq11 by default)",
        type=str, default=None)

    parser.add_argument(
        "--filename", help="an output filename", type=str, default=None)

    req_args.add_argument(
        "--out", help="path to output directory", type=str, required=True)

    return parser


# -----------------------------------------------------------------------------
if __name__ == '__main__':

    # call usage and parse arguments
    arg_parser = usage()
    args = arg_parser.parse_args()

    seq_file = ''
    filename = ''
    sql = ''

    # set sequence file to point to rfamseq
    if args.infile is None:
        seq_file = rfam_config.RFAMSEQ_PATH

    # check valid destination directory
    if not os.path.isdir(args.out):
        print "\nIncorrect output directory. "
        sys.exit()

    if args.filename is None:
        filename = OUT_FILE_NAME
    else:
        filename = args.filename

    export_sequences(seq_file, sql, filename, args.out)
