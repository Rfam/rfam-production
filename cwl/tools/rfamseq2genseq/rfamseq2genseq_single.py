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

# -----------------------------------------------------------------------------

import os
import sys


# -----------------------------------------------------------------------------


def convert_rfamseq_to_genseq(rfamseq_file):
    """
    Converts an rfamseq file to genseq to map genome (upid) and sequence
    accessions

    :param rfamseq_file: A genome specific rfamseq file in the form of
    upid.rfamseq, as generated from rfamseq table

    returns: void
    """
    # store output in input file directory
    dest_dir = os.path.split(rfamseq_file)[0]
    # get the input filename to keep track of the upids
    filename = os.path.basename(rfamseq_file).partition('.')[0]
    genseq_file = open(filename + '.genseq', 'w')
    rfamseq_fp = open(rfamseq_file, 'r')
    for line in rfamseq_fp:
        line = line.strip().split('\t')
        genseq_file.write(filename + '\t' + line[0] + '\n')
    genseq_file.close()
    rfamseq_fp.close()


if __name__ == '__main__':
    rfamseq_file = sys.argv[1]
    convert_rfamseq_to_genseq(rfamseq_file)
