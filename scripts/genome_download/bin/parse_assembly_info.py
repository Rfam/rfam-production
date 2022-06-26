#!/usr/bin/env python3

import csv

from sqlitedict import SqliteDict
import click

@click.command()
@click.argument('info_file', type=click.File('r'))
@click.argument('output', type=click.Path())
def main(info_file, output):
    reader = csv.DictReader(info_file, delimiter='\t')
    summary = {}
    with SqliteDict(output) as db:
        for index, info in enumerate(reader):
            accession = info['assembly_accession']
            prev = summary.get(accession, None)
            if prev is not None and prev != info:
                raise ValueError("Mutiple mappings to %s" % accession)
            db[accession] = info
            if index % 100000 == 0:
                db.commit()
    db.commit()

if __name__ == '__main__':
    main()
