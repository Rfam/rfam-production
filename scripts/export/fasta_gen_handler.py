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
Description:    Calls fasta_generator to generate fasta files for all Rfam
                families in rfam_live

Comments:       It is a prerequisite that the sequence file is indexed using
                esl-sfetch --index option
"""

# ---------------------------------IMPORTS-------------------------------------

import os
import sys
import subprocess
from utils import RfamDB
from config import rfam_config

# -----------------------------------------------------------------------------

LSF_GROUP = rfam_config.FA_EXPORT_GROUP

# -----------------------------------------------------------------------------


def fasta_gen_handler(seq_file, out_dir):
    """
    The purpose of this script is to handle the fasta generation process,
    generate individual shell scripts for each available family and submit
    them to the cluster

    seq_file:   Path to the input sequence file (e.g. rfamseq11.fa)
    out_dir:    The output directory where the fasta files will be generated

    """

    # fetch family accessions
    cnx = RfamDB.connect()

    cursor = cnx.cursor(buffered=True)

    query = ("SELECT rfam_acc FROM family")

    cursor.execute(query)

    families = cursor.fetchall()

    cursor.close()
    RfamDB.disconnect(cnx)

    # create scripts dir within output directory
    if (not os.path.exists(os.path.join(out_dir, "scripts"))):
        os.mkdir(os.path.join(out_dir, "scripts"))

    if (not os.path.exists(os.path.join(out_dir, "log"))):
        os.mkdir(os.path.join(out_dir, "log"))

    for fam in families:

        # 1. Generate script file
        sh_path = shell_script_generator(
            seq_file, str(fam[0]), out_dir, os.path.join(out_dir, "scripts"))

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

    fp = open(file_path, 'w')

    fp.write("#!/bin/csh\n")
    fp.write("#BSUB -q research-rh7\n")
    fp.write("#BSUB -M 8000\n")
    fp.write("#BSUB -R \"rusage[mem=8000,tmp=1000]\"\n")
    fp.write("#BSUB -o \"/tmp/%J.out\"\n")
    fp.write("#BSUB -e \"/tmp/%J.err\"\n")

    fp.write(
        "#BSUB -f \"%s/%s.out < /tmp/%sJ.out\"\n" % (log_dir, rfam_acc, chr(37)))

    fp.write(
        "#BSUB -f \"%s/%s.err < /tmp/%sJ.err\"\n" % (log_dir, rfam_acc, chr(37)))

    fp.write("#BSUB -Ep \"rm /tmp/$LSB_JOBID.*\"\n")
    fp.write("#BSUB -g %s \n\n" % (LSF_GROUP))
    fp.write("python %s %s %s %s \n" %
             (rfam_config.FA_GEN, seq_file, rfam_acc, fa_outdir))

    fp.close()

    return file_path

# -----------------------------------------------------------------------------


def usage():
    """
    Displays information on how to run fasta_gen_handler
    """

    print "\nUsage:\n------"

    print "\npython fasta_gen_handler.py seq_file out_dir"

    print "\nseq_file: Path to sequence for sequence export (e.g. rfamseq11.fa)"
    print "out_dir: The path to the output directory"
    print "\n-h option to display usage\n"

# -----------------------------------------------------------------------------

if __name__ == '__main__':

    # minor input checks
    if (sys.argv[1] == "-h"):
        usage()
        sys.exit()

    elif(len(sys.argv) == 3):
        seq_file = sys.argv[1]
        out_dir = sys.argv[2]

        if (os.path.isfile(seq_file) and os.path.isdir(out_dir)):
            fasta_gen_handler(seq_file, out_dir)

        else:
            print "\nIncorrect Input."
            usage()

    else:
        usage()
