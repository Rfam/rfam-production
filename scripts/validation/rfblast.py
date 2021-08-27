#!/usr/bin/env python

"""
Replace unknown accessions in SEED alignment using NCBI BLAST results.

Usage:
cd /path/to/folder/with/SEED
rfblast.py XXXXX

where XXXXX.json is the NCBI BLAST result.
"""


import json
import os

import click

IDENTITY = 90
QUERY_COVERAGE = 70

#TODO do not replace valid accessions
#TODO create tests using the H742DVSC01R-Alignment.json file

def get_accession(gid):
    """
    Get versioned accession, for example:
    Input:
    gi|2047803076|gb|CP061286.1|
    Output:
    CP061286.1
    """
    parts = gid.split('|')
    return parts[-2]


def get_blast_data(filename):
    """
    Load BLAST JSON output.
    """
    with open(filename, 'r') as f:
        return json.load(f)


def choose_replacement(data, min_identity, min_query_coverage):
    """
    Loop over BLAST results and pick best replacement for each hit.
    """
    # do not pick replacement from the same accession if already seen
    fasta = ''
    seen_accessions = set()
    for query_num, search in enumerate(data['BlastOutput2']):
        query_title = search['report']['results']['search']['query_title']
        query_len = search['report']['results']['search']['query_len']
        replacement_found = False
        for entry in search['report']['results']['search']['hits']:
            acc = get_accession(entry['description'][0]['id'])
            if acc not in seen_accessions:
                seen_accessions.add(acc)
                replacement_found = True
            else:
                continue
            sequence = entry['hsps'][0]['hseq']
            start = entry['hsps'][0]['hit_from']
            end = entry['hsps'][0]['hit_to']
            align_len = entry['hsps'][0]['align_len']
            gaps = entry['hsps'][0]['gaps']
            exact_matches = entry['hsps'][0]['identity']
            identity = float(exact_matches) / align_len * 100
            query_coverage = float(align_len - gaps) / query_len * 100
            target_coverage = float(align_len - gaps) / len(sequence) * 100
            if identity >= min_identity and query_coverage >= min_query_coverage:
                warning = ''
            else:
                warning = '      WARNING: '
            summary = ('#{query_num} {warning}Replace {query_title} '
                       'with {acc}/{start}-{end} at {identity}% identity; '
                       '{gaps} gaps; query coverage {query_coverage}').format(
                acc=acc, start=start, end=end, query_title=query_title,
                identity=round(identity), query_coverage=round(query_coverage),
                target_coverage=round(target_coverage, 2), gaps=gaps,
                warning=warning,
                query_num=query_num+1
            )
            print(summary)
            if not warning:
                fasta += '>{acc}/{start}-{end}\n{sequence}\n'.format(
                    acc=acc,
                    start=start,
                    end=end,
                    sequence=sequence.replace('-', '').replace('T', 'U')
                )
            if replacement_found:
                break
    return fasta


def generate_new_seed(fasta, destination):
    filename = 'replacement.fasta'
    replacement_fasta = os.path.join(destination, filename)
    with open(replacement_fasta, 'w') as f:
        f.write(fasta)
    cmd = ('cd {destination} && '
           'cmalign --noprob CM {filename} > tempseed && '
           'esl-reformat pfam tempseed > NEWSEED && '
           'echo "Old SEED info:" && esl-alistat SEED && '
           'echo "New SEED info:" && esl-alistat NEWSEED && '
           'rm tempseed && cd -').format(destination=destination, filename=filename)
    os.system(cmd)


@click.command()
@click.option('--identity', default=IDENTITY, help='Minimum % identity between query and target')
@click.option('--query_coverage', default=QUERY_COVERAGE, help='Minimum coverage of the seed sequence')
def rfblast(identity, query_coverage):
    filename = 'carA-HHBXHPJ5013-Alignment.json'
    destination = 'temp/carA'
    blast_data = get_blast_data(filename)
    fasta = choose_replacement(blast_data, identity, query_coverage)
    generate_new_seed(fasta, destination)


if __name__ == '__main__':
    rfblast()
