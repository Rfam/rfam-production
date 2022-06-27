#!/usr/bin/env python3

import csv

from sqlitedict import SqliteDict
import click

@click.command()
@click.argument('info_file', type=click.File('r'))
@click.argument('output', type=click.Path())
def main(info_file, output):
    count = 0
    reader = csv.DictReader(info_file, delimiter='\t')
    with SqliteDict(output) as db:
        for index, info in enumerate(reader):
            count += 1
            accession = info['assembly_accession']
            if accession in db:
                raise ValueError("Mutiple mappings to %s" % accession)
            db[accession] = info
            if index % 100000 == 0:
                db.commit()
        db.commit()

    with SqliteDict(output) as db:
        assert len(db) == count, "Did not load all assembiles"

if __name__ == '__main__':
    main()
