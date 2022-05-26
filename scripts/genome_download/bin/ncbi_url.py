#!/usr/bin/env python

import pickle

import click

@click.command()
@click.argument('ncbi-file', type=click.File('rb'))
@click.argument('gca')
def main(ncbi_file, gca):
    ncbi = pickle.load(ncbi_file)
    if gca not in ncbi:
        raise ValueError("Unknown GCA id %s" % gca)
    info = ncbi[gca]
    path = info['ftp_path']
    parts = path.split('/')
    name = parts[-1]
    url = f"{path}/{name}_genomic.fna.gz"
    print(url)


if __name__ == '__main__':
    main()
