#!/usr/bin/env python3

import csv
import pickle

import click

@click.command()
@click.argument('info_file', type=click.File('r'))
@click.argument('output', type=click.File('wb'))
def main(info_file, output):
    reader = csv.DictReader(info_file, delimiter='\t')
    summary = {}
    for info in reader:
        accession = info['assembly_accession']
        prev = summary.get(accession, None)
        if prev is not None and prev != info:
            raise ValueError("Mutiple mappings to %s" % accession)
        summary[accession] = info
    pickle.dump(summary, output)

if __name__ == '__main__':
    main()
