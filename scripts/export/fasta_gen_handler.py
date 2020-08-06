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
Description:    Calls fasta_generator to generate fasta files for all Rfam
                families in rfam_live

Comments:       It is a prerequisite that the sequence file is indexed using
                esl-sfetch --index option
"""

# ---------------------------------IMPORTS-------------------------------------

import os
import sys
import subprocess
import argparse
from utils import RfamDB
from config import rfam_config

# -----------------------------------------------------------------------------

LSF_GROUP = rfam_config.FA_EXPORT_GROUP

# -----------------------------------------------------------------------------


def fasta_gen_handler(seq_file, out_dir, rfam_accessions=None):
    """
    The purpose of this script is to handle the fasta generation process,
    generate individual shell scripts for each available family and submit
    them to the cluster

    seq_file:   Path to the input sequence file (e.g. rfamseq11.fa)
    out_dir:    The output directory where the fasta files will be generated

    """

    # fetch family accessions
    families = []

    if rfam_accessions is None:
        cnx = RfamDB.connect()

        cursor = cnx.cursor(buffered=True)

        query = ("SELECT rfam_acc FROM family")

        cursor.execute(query)

        entries = cursor.fetchall()

        cursor.close()
        RfamDB.disconnect(cnx)

        families = [str(fam[0]) for fam in entries]
    else:
        fp = open(rfam_accessions, 'r')
        families = [x.strip() for x in fp]
        fp.close()

    # create scripts dir within output directory
    if not os.path.exists(os.path.join(out_dir, "scripts")):
        os.mkdir(os.path.join(out_dir, "scripts"))

    if not os.path.exists(os.path.join(out_dir, "log")):
        os.mkdir(os.path.join(out_dir, "log"))

    for fam in families:

        # 1. Generate script file
        sh_path = shell_script_generator(
            seq_file, fam, out_dir, os.path.join(out_dir, "scripts"))

        # 2. submit job under group
        cmd = "bsub < %s" % (sh_path)
        subprocess.call(cmd, shell=True)

# -----------------------------------------------------------------------------


def shell_script_generator(seq_file, rfam_acc, fa_outdir, out_dir=None):
    """
    Generates family specific shell scripts to split fasta generation into
    individual jobs

    seq_file:   The path to sequence file (e.g. )
    rfam_acc:   A valid Rfam family accession
    fa_outdir:  A path to where fasta files will be generated
    out_dir:    A path to an output directory where the shell scripts will be
                generated. If None, fa_outdir is used by default
    """

    # If no specific directory is provided for the shell scripts, generate them
    # in the fa output directory

    file_path = ''
    if out_dir is None:
        file_path = os.path.join(fa_outdir, rfam_acc + ".sh")

    else:
        file_path = os.path.join(out_dir, rfam_acc + ".sh")

    log_dir = os.path.join(fa_outdir, "log")

    output_fp = open(file_path, 'w')

    output_fp.write("#!/bin/csh\n")
    output_fp.write("#BSUB -q research-rh74\n")
    output_fp.write("#BSUB -M 8000\n")
    output_fp.write("#BSUB -R \"rusage[mem=8000,tmp=1000]\"\n")
    output_fp.write("#BSUB -o \"/tmp/%J.out\"\n")
    output_fp.write("#BSUB -e \"/tmp/%J.err\"\n")

    output_fp.write(
        "#BSUB -f \"%s/%s.out < /tmp/%sJ.out\"\n" % (log_dir, rfam_acc, chr(37)))

    output_fp.write(
        "#BSUB -f \"%s/%s.err < /tmp/%sJ.err\"\n" % (log_dir, rfam_acc, chr(37)))

    output_fp.write("#BSUB -Ep \"rm /tmp/$LSB_JOBID.*\"\n")
    output_fp.write("#BSUB -g %s \n\n" % (LSF_GROUP))
    output_fp.write("python %s %s %s %s \n" % (rfam_config.FA_GEN, seq_file,
                                               rfam_acc, fa_outdir))

    output_fp.close()

    return file_path

# -----------------------------------------------------------------------------


def parse_arguments():
    """
    Basic argument parsing using python's argparse

    return: Argparse parser object
    """

    parser = argparse.ArgumentParser("Rfam fasta file generation handler")

    parser.add_argument('--seq-db', help="Sequence database in fasta format",
                        action="store", default=None)
    parser.add_argument('--outdir', help="Output directory", action="store")

    parser.add_argument('-f', help='A file with a list of Rfam family accessions',
                        action="store", default=None)

    return parser

# -----------------------------------------------------------------------------

if __name__ == '__main__':

    parser = parse_arguments()
    args = parser.parse_args()

    seq_db = args.seq_db
    outdir = args.outdir

    if os.path.isfile(seq_db) and os.path.isdir(outdir):
        fasta_gen_handler(seq_db, outdir, args.f)


