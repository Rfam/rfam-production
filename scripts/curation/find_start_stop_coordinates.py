#!/usr/bin/env python
"""
Find start/stop coordinates given a Stockholm file with valid
INSDC or RefSeq accessions.

Before:
AE015928 GUGUUUUUCAUAGUA
AP006841 GUGUUUUUCAUAGUA
CP000139 GUGUUUUUCAUAGUA

After:
AE015928.1/959075-959219     GUGUUUUUCAUAGUA
AP006841.1/2613625-2613766   GUGUUUUUCAUAGUA
CP000139.1/3462469-3462596   GUGUUUUUCAUAGUA

Note that the input file needs to contain valid accessions with or without
versions. The accessions cannot have any other suffixes - these need to be
removed manually before running the script.

Usage:
python find_start_stop_coordinates.py <input.sto>

The command generates a new file <input_with_accessions.sto>
"""

import os
import re
import sys
import time

from Bio import AlignIO
from Bio import SeqIO


def get_accession_version(filename):
    """
    Get version of an accession based on a fasta file retrieved from NCBI.

    For example, in the following line the version is 2
    >NC_034265.1 Tobacco virus 2, complete genome
    """
    version = '.1'
    with open(filename, 'r') as f_in:
        for line in f_in:
            if not line.startswith('>'):
                print('Invalid fasta file')
                return
            fields = line.split(' ')
            name = fields[0]
            _, version = name.split('.')
            version = '.' + version
            break
    return version


def get_fasta_file(identifier):
    """
    Download fasta file from NCBI if not already downloaded
    and store in a folder.
    """
    script_location = os.path.dirname(os.path.abspath(__file__))
    fasta_folder = os.path.join(script_location, 'fasta')
    os.system('mkdir -p {}'.format(fasta_folder))
    fasta_file = os.path.join(fasta_folder, identifier + '.fasta')
    if not os.path.exists(fasta_file) or os.stat(fasta_file).st_size == 0:
        time.sleep(1)
        cmd = 'wget -O {0} "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={1}&rettype=fasta&retmode=text"'.format(fasta_file, identifier)
        os.system(cmd)
    return fasta_file


def map_accessions(filename):
    print(filename)
    mapping = {}
    align = AlignIO.read(filename, 'stockholm')
    for record in align:
        fasta_file = get_fasta_file(record.id)
        fasta = SeqIO.read(fasta_file,"fasta")
        if '.' not in record.id:
            version = get_accession_version(fasta_file)
        else:
            version = ''

        s1 = record.seq.ungap('.').ungap('-').upper().back_transcribe()
        start = fasta.seq.find(s1)
        if start == -1:  # not found in forward orientation
            s2 = s1.reverse_complement()  # try reverse orientation
            start = fasta.seq.find(s2)
            if start != -1:
                new_accession = '{}{}/{}-{}'.format(record.id, version, start + len(s1), start + 1)
        else:
            new_accession = '{}{}/{}-{}'.format(record.id, version, start + 1, start + len(s1))

        if start != -1:
            print('{} corresponds to {}'.format(record.id, new_accession))
            mapping[record.id] = new_accession
        else:
            print('\nWarning: sequence {} not found and will be removed\n'.format(record.id))
    return mapping


def write_new_seed(in_filename, out_filename, id_mapping):
    """
    Biopython ignores some metadata so cannot be used to output
    an updated version of the alignment.
    """
    old_accs = id_mapping.keys()
    with open(in_filename, 'r') as f_in:
        with open(out_filename, 'w') as f_out:
            for line in f_in:
                if line.startswith('# STOCKHOLM'):
                    f_out.write(line)
                if line.startswith('//') or line == "\n":
                    f_out.write(line)
                if line.startswith('#'): #GR or #GS lines
                    parts = re.split(r'\s+', line.strip())
                    info = ' '.join(parts[:-1])
                    fixed_width_info = '{0: <30}'.format(info)
                    f_out.write(fixed_width_info + parts[-1] + '\n')
                for old_acc in old_accs:
                    if old_acc not in id_mapping:
                        print('Skipping {}'.format(old_acc))
                        continue
                    if line.startswith(old_acc):
                        fixed_width_acc = '{0: <30}'.format(id_mapping[old_acc])
                        parts = re.split(r'\s+', line.strip())
                        f_out.write(fixed_width_acc + parts[-1] + '\n')
    print('Generated new file {}'.format(out_filename))


def get_new_file_name(in_filename):
    """
    Given a path to an input file, generate a path to a new file
    that ends with "_with_accessions".
    """
    head, tail = os.path.split(in_filename)
    extension = os.path.splitext(in_filename)[1]
    new_name = tail.replace(extension, '_with_accessions' + extension)
    return os.path.join(head, new_name)


def main():
    if len(sys.argv) != 2:
        raise Exception('Provide a Stockholm file as an argument')
    in_filename = sys.argv[1]
    out_filename = get_new_file_name(in_filename)
    id_mapping = map_accessions(in_filename)
    write_new_seed(in_filename, out_filename, id_mapping)
    cmd = 'esl-alistat {}'.format(out_filename)
    os.system(cmd)


if __name__ == '__main__':
    main()
