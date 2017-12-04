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
Usage:

python automate_rfmake.py -i /path/to/emerge_file -d /path/to/output
"""


import argparse
import os
import sys


def parse_input_file(filename):
    """
    Read input data and standardise sequences and names.
    Example input is provided in `example.tsv`.
    """
    with open(filename, 'r') as f:
        lines  = f.readlines()
        return [x.strip() for x in lines]

def run(args):
    """
    * run rfmake
    * check for overlaps
    """
    for rna in parse_input_file(args.inputfile):
        folder_name = os.path.join(args.destination, '%s' % rna)
        if not os.path.exists(folder_name):
            os.mkdir(folder_name)
        os.chdir(folder_name)
        cmd = ('bsub -o {0}/lsf_rfmake_output.txt -e {0}/lsf_rfmake_error.txt -g /emerge '
                     '"cd {0} && '
                     'rfmake.pl -t 50 -a && cd .. && '
                     'rqc-overlap.pl {0}"').format(folder_name)
        print cmd
        os.system(cmd)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--destination', default=os.getcwd(), help='Specify folder where the output will be created')
    parser.add_argument('-i', '--inputfile', default='example.tsv', help='Specify input file with names of folders where the analysis needs to be done')
    args = parser.parse_args()

    if not args.inputfile:
        print 'Please specify input file'
        sys.exit()

    run(args)
