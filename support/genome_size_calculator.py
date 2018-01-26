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

# --------------------------------------------------------------------------------------------------

import os
import sys
import subprocess
import json
from config import gen_config as gc

# --------------------------------------------------------------------------------------------------


def calculate_genome_size(genome_source):
    """
    Calculate the size of a given sequence file (this should be a genome) and return the total
    number of nucleotides in the file, or if genome_fasta is a directory, return a dictionary
    of genome_ids: size pairs

    genome_source: A valid sequence file in fasta format or a directory of fasta files

    return: A dictionary or the an integer number which is the total number of nucleotides
    found in the sequence file that is provided as input
    """

    genome_files = []
    genome_sizes = {}

    if os.path.isdir(genome_source):
        genome_files = [x for x in os.listdir(genome_source)
                        if x.endswith('.fa') or x.endswith('.fasta')]

        for seq_file in genome_files:
            try:
                args = [gc.ESL_SEQSTAT, '--dna', os.path.join(genome_source, seq_file)]
                # response = subprocess.check_output(args, shell=True, stderr=subprocess.STDOUT)
                # calculate genome size using esl-seqstat and return Total
                response = subprocess.check_output(args)
                total_res_line = response.split('\n')[3]

            except subprocess.CalledProcessError as err:
                raise RuntimeError("command '{}' return with error (code {}): {}".format(err.cmd,
                                                                                         err.returncode,
                                                                                         err.output))

            # split the line and get total number of residues which must be the last element in the list
            # then strip any white characters and convert to integer
            total_res = int(total_res_line.split(' ')[-1].strip())
            genome_id = seq_file.partition('.')[0]
            genome_sizes[genome_id] = total_res

    else:

        try:
            # response = subprocess.check_output(args, shell=True, stderr=subprocess.STDOUT)
            # calculate genome size using esl-seqstat and return Total
            args = [gc.ESL_SEQSTAT, '--dna', genome_source]
            response = subprocess.check_output(args)
            total_res_line = response.split('\n')[3]

        except subprocess.CalledProcessError as err:
            raise RuntimeError("command '{}' return with error (code {}): {}".format(err.cmd,
                                                                                     err.returncode,
                                                                                     err.output))

        # split the line and get total number of residues which must be the last element in the list
        # then strip any white characters and convert to integer
        total_res = int(total_res_line.split(' ')[-1].strip())
        return total_res

    return genome_sizes

# --------------------------------------------------------------------------------------------------


if __name__ == '__main__':

    # first argument is either a sequence file or a directory of multiple sequence files
    # second argument is the name of the output filename

    if len(sys.argv) == 3:
        if os.path.isdir(sys.argv[1]):
            output_fp = open(os.path.join(sys.argv[1], sys.argv[2]+'.json'), 'w')
            genome_sizes = calculate_genome_size(sys.argv[1])
            json.dump(genome_sizes, output_fp)
            output_fp.close()

    # print genome size as total # of nucleotides
    else:
        print calculate_genome_size(sys.argv[1])


