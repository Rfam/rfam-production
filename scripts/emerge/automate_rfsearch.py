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

python automate_rfsearch.py -i /path/to/emerge_file -d /path/to/output

Analysing result:

# check that all LSF jobs completed successfully
ls | wc -l
find . -type f -name lsf_output.txt | xargs grep 'Success' | wc -l

# count the number of lines above the best reversed hit
find . -type f -name outlist -exec sh -c 'sed -n "0,/REVERSED/p" {} | wc -l' \; -print
"""


import argparse
import csv
import os
import sys


def parse_input_file(filename):
    """
    Read input data and standardise sequences and names.
    Example input is provided in `example.tsv`.
    """
    with open(filename, 'r') as tsv:
        reader = csv.DictReader(tsv, delimiter='\t')
        for row in reader:
            sequence = row['Sequence (RNA or DNA)'].replace('-', '').upper()
            if len(sequence) < 50:
                continue
            yield {
                'sequence': sequence,
                'name': row['ncRNA name'].replace(' ', '_').replace("'", ''),
                'row_id': row['No.'],
            }

def run(args):
    """
    * create FASTA file
    * predict secondary structure
    * make SEED
    * launch rfsearch
    """
    for rna in parse_input_file(args.inputfile):
        folder_name = os.path.join(args.destination, '%s_%s' % (rna['row_id'], rna['name']))
        if not os.path.exists(folder_name):
            os.mkdir(folder_name)
        else:
            outlist = os.path.join(folder_name, 'outlist')
            if os.path.exists(outlist):
                continue
        os.chdir(folder_name)
        with open('input.fasta', 'w') as fasta:
            fasta.write('>%s\n%s\n' % (rna['name'], rna['sequence']))
        cmd = ('bsub -o {0}/lsf_output.txt -e {0}/lsf_error.txt -g /emerge '
                     '"cd {0} && '
                     'predict_ss.pl -infile input.fasta -outfile SEED -r && '
                     'rfsearch.pl -nodesc -t 30 -cnompi -relax"').format(folder_name)
        print cmd
        os.system(cmd)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--destination', default=os.getcwd(), help='Specify folder where the output will be created')
    parser.add_argument('-i', '--inputfile', default='example.tsv', help='Specify input file with E-Merge data')
    args = parser.parse_args()

    if not args.inputfile:
        print 'Please specify input file'
        sys.exit()

    run(args)
