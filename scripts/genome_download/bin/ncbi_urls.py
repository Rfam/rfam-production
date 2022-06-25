#!/usr/bin/env python

import json
import pickle
from pathlib import Path

import click


@click.command()
@click.argument('ncbi-file', type=click.File('rb'))
@click.argument('gca_file', type=click.File('r'))
@click.argument('directory', type=click.Path())
def main(ncbi_file, gca_file, directory):
    ncbi = pickle.load(ncbi_file)
    base = Path(directory)
    for line in gca_file:
        row = json.loads(line)
        gca = row['accession']
        if gca not in ncbi:
            raise ValueError("Unknown GCA id %s" % gca)
        info = ncbi[gca]
        path = info['ftp_path']
        parts = path.split('/')
        name = parts[-1]
        save_path = base / f"{row['upi']}.fa.gz"
        url = f"{path}/{name}_genomic.fna.gz"
        print(save_path)
        print(url)


if __name__ == '__main__':
    main()
