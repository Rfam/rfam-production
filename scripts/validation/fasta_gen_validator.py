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
Description:    Validation script to check fasta file generation process
"""

# ---------------------------------IMPORTS-------------------------------------

import os
import sys
import gzip
import string
from utils import RfamDB

# -----------------------------------------------------------------------------


def get_full_region_seq_counts():
    """
    Builds a dictionary where keys are Rfam family accessions (rfam_acc)
    and values the number of sequences in full_region per family
    (e.g. {'RFXXXXX':N,...})
    """

    seq_counts = {}

    # get a connection object for RfamDB
    cnx = RfamDB.connect()

    cursor = cnx.cursor(buffered=True)

    query = ("SELECT rfam_acc, count(*) FROM full_region\n"
             "GROUP BY rfam_acc")

    cursor.execute(query)

    # get full_region sequence counts per family
    raw_counts = cursor.fetchall()

    # build dictionary
    for entry in raw_counts:
        seq_counts[str(entry[0])] = int(entry[1])

    # close DB handles
    cursor.close()
    RfamDB.disconnect(cnx)

    # result dictionary
    return seq_counts

# -----------------------------------------------------------------------------


def get_fasta_seq_counts(fasta_files):
    """
    This module uses database data to check whether the fasta generation
    process was successful and returns a list of the families that need to
    be processed again

    fasta_files: The path to the fasta files directory
    """

    # work on fasta files
    fa_files = os.listdir(fasta_files)

    fa_files = filter(
        lambda x: x.endswith(".fa") or x.endswith(".fa.gz"), fa_files)

    fa_seq_counts = {}
    count = 0

    # get sequence counts per family
    for fa in fa_files:
        # get rfam_acc
        rfam_acc = string.split(fa, '.')[0]
        fp = gzip.open(os.path.join(fasta_files, fa), 'r')

        for line in fp:
            if line[0] == '>':
                count = count + 1

        # update dictionary
        fa_seq_counts[rfam_acc] = count

        count = 0
        fp.close()

    return fa_seq_counts

# -----------------------------------------------------------------------------


def compare_seq_counts(db_counts, fa_counts):
    """
    Compares the number of sequences per family in full_region table with
    the number of sequences written in the distinct fasta files

    db_counts:  A dictionary with the number of sequences per family as
                found in full_region (e.g. {'RFXXXXX':N,...}). Output of
                get_full_region_seq_counts
    fa_counts:  A dictionary with the number of sequences per family fasta
                file (e.g. {'RFXXXXX':N,...}). Output of get_fasta_seq_counts
    """

    faulty_fams = []

    for rfam_acc in db_counts.keys():

        if (db_counts[rfam_acc] != fa_counts[rfam_acc]):
            faulty_fams.append(rfam_acc)

    return faulty_fams

# -----------------------------------------------------------------------------


def usage():
    """
    Displays information on how to run fasta_gen_validator
    """
    print "\nUsage:\n------"
    print "\npython fasta_gen_validator.py /path/to/fasta_files"
    print "\nfasta_files: The path to the fasta files directory\n"

# -----------------------------------------------------------------------------

if __name__ == '__main__':

    # path to fasta files directory
    fasta_files = sys.argv[1]

    # some input checks
    if (not os.path.isdir(fasta_files)):
        print "\nIncorrect Input."
        usage()
        sys.exit()

    elif(sys.argv[1] == "-h"):
        usage()
        sys.exit()

    else:
        # get sequence counts from DB
        db_seq_counts = get_full_region_seq_counts()

        # get sequence counts from fasta_files
        fa_seq_counts = get_fasta_seq_counts(fasta_files)

        # compare counts and print families that need to be re-submitted to lsf
        faulty_fams = compare_seq_counts(db_seq_counts, fa_seq_counts)

        # export completed successfully
        if (len(faulty_fams) == 0):
            print "\nFasta export process completed successfully!"

        # print family accessions with erroneous files
        else:
            print "\nErroneous files: "
            for rfam_acc in faulty_fams:
                print rfam_acc
