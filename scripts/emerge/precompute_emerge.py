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

1. Make sure that rfsearch, rfmake and other commands are in PATH.

2. python automate_rfsearch.py -i /path/to/emerge_file -d /path/to/output

Analysing result:

# check that all LSF jobs completed successfully
ls | wc -l
find . -type f -name lsf_output.txt | xargs grep 'Success' | wc -l

# find jobs that didn't finish successfully
find . -type f -name lsf_output.txt | xargs grep -L 'Success'

# find all directories without the outlist file
find . -maxdepth 1 -mindepth 1 -type d | while read dir; do [[ ! -f $dir/outlist ]] && echo "$dir has no outlist"; done

# count the number of lines above the best reversed hit
find . -type f -name outlist -exec sh -c 'sed -n "0,/REVERSED/p" {} | wc -l' \; -print

# get overlapping Rfam families
find . -type f -name overlap | xargs grep -o -P "RF\d{5}" | sort | uniq
"""


import argparse
import csv
import os
import re
import sys


def parse_input_file(filename):
    """
    Read input data and standardise sequences and names.
    Example input is provided in `example.tsv`.
    """
    skipped = {
        'length': 0,
        'in_rfam': 0,
    }
    MIN_LENGTH = 50
    with open(filename, 'r') as tsv:
        reader = csv.DictReader(tsv, delimiter='\t')
        for row in reader:
            sequence = row['Sequence (RNA or DNA)'].replace('-', '').replace('.','').upper()
            sequence = re.sub(r'\s', '', sequence)
            if len(sequence) < MIN_LENGTH:
                skipped['length'] += 1
                continue
            if row['In Rfam? http://rfam.xfam.org/search'].startswith('RF0'):
                skipped['in_rfam'] += 1
                continue
            yield {
                'sequence': sequence,
                'name': row['ncRNA name'].replace(' ', '_').replace("'", '').replace('/', '-').replace('(', '-').replace(')', '-'),
                'row_id': row['No.'],
            }
        print 'Skipped %i sequences shorter than %i nucleotides' % (skipped['length'], MIN_LENGTH)
        print 'Skipped %i sequences already in Rfam' % (skipped['in_rfam'])

def run(args):
    """
    * create FASTA file
    * predict secondary structure
    * make SEED
    * launch rfsearch
    """
    for rna in parse_input_file(args.inputfile):
        folder = '%s_%s' % (rna['row_id'], rna['name'])
        rna_dir = os.path.join(args.destination, folder)
        if not os.path.exists(rna_dir):
            os.mkdir(rna_dir)
        else:
            overlap = os.path.join(rna_dir, 'overlap')
            if os.path.exists(overlap):
                continue
        os.chdir(rna_dir)
        with open('input.fasta', 'w') as fasta:
            fasta.write('>%s\n%s\n' % (rna['name'], rna['sequence']))
        cmd = ('module load mpi/openmpi-x86_64 && '
               'bsub -o {0}/lsf_output.txt -e {0}/lsf_error.txt -g /emerge '
                     '"cd {0} && '
                     'rm -f DESC && '
                     'predict_ss.pl -infile input.fasta -outfile SEED -r && '
                     'rfsearch.pl -nodesc -t 30 -relax && '
                     'rfmake.pl -t 50 -a -forcethr && '
                     'mkdir rscape && R-scape --outdir rscape --cyk align && '
                     'cd .. && '
                     'rqc-overlap.pl {1}"').format(rna_dir, folder)
        print cmd
        if not args.test:
            os.system(cmd)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--destination', default=os.getcwd(), help='Specify folder where the output will be created')
    parser.add_argument('-i', '--inputfile', default='example.tsv', help='Specify input file with E-Merge data')
    parser.add_argument('-t', '--test', action='store_true', help='Test mode: print commands and exit')
    parser.set_defaults(test=False)
    args = parser.parse_args()

    if not args.inputfile:
        print 'Please specify input file'
        sys.exit()

    run(args)
