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
import os
import re
import subprocess
from subprocess import PIPE
from config import rfam_local as rl

# --------------------------------------------------------------------------------------------------

def seq_validator(sequence):
    """
    Checks if the sequence provided is valid fasta sequence. Returns True
    if the sequence is valid, otherwise returns False.

    sequence: A string for validation
    """

    # checks for ascii characters that should not appear in a fasta sequence
    seq_val = re.compile(r"[.-@|\s| -)|z-~|Z-`|EFIJLOPQX|efijlopqx+,]+")

    if(seq_val.search(sequence) is None):
        return True

    return False

# --------------------------------------------------------------------------------------------------
def get_nt_count(seq_file, type="dna"):
    """
    Using esl-seqstat

    type: A string any of "dna", "rna", "amino"
    :return:
    """
    nt_count = 0
    inf_type = "--%s"%type
    p = subprocess.Popen([rl.ESL_SEQSTAT, inf_type, seq_file], stdout=PIPE, stderr=PIPE)
    data = p.communicate()[0]

    values = data.split('\n')
    for stat in values:
        if stat.find("Total") != -1:
            stat = stat.split(' ')
            nt_count = int(stat[len(stat) - 1])

    return nt_count

# --------------------------------------------------------------------------------------------------
def calculate_seqdb_size(seqdb, type="dna"):
    """
    This function uses infernal's easel package to
    :param seqdb: fasta file or a directory of fasta files
    :return:
    """
    db_size = 0
    #single seqdb file
    if os.path.isfile(seqdb):
        db_size = get_nt_count(seqdb, type)

    #directory
    else:
        seq_files = [x for x in os.listdir(seqdb) if x.endswith(".fa")]
        for seq_file in seq_files:
            seq_file_size = get_nt_count(os.path.join(seqdb, seq_file), type)
            db_size += seq_file_size

    return db_size

# --------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    #seqdb = "/Users/ikalvari/Desktop/tests/gen_dwl/euk/UP000005640/KI270724.1.fa"
    #seqdb = "/Users/ikalvari/Desktop/tests/gen_dwl/euk/UP000005640/CM000672.2.fa"
    seqdb = "/Users/ikalvari/Desktop/tests/gen_dwl/euk/UP000005640"


    nt_count = calculate_seqdb_size(seqdb, "dna")
    #nt_count = nucleotide_counter(seqdb)


    print nt_count

