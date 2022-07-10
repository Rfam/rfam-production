#!/usr/bin/env python

import re
import json
from pathlib import Path
import typing as ty

import click
from sqlitedict import SqliteDict


def ftp_path(info: SqliteDict, accession: str) -> ty.Optional[str]:
    path = None
    if accession not in info:
        return None
    else:
        path = info[accession]['ftp_path']

    parts = path.split('/')
    name = parts[-1]
    return f"{path}/{name}_genomic.fna.gz"


def determine_latest(db: SqliteDict, id: str) -> str:
    if id in db:
        return id
    possible = {}
    pattern = re.compile(f"^{id}.(\d+)$")
    for key in db.iterkeys():
        if match := re.match(pattern, key):
            index = int(match.group(1))
            possible[index] = key
    to_use = max(possible.keys())
    return possible[to_use]


@click.command()
@click.argument('ncbi-file', type=click.Path())
@click.argument('gca_file', type=click.File('r'))
@click.argument('directory', type=click.Path())
@click.argument('url-file', default='-', type=click.File('w'))
@click.argument('ena-only', default='ena-only.jsonl', type=click.File('w'))
def main(ncbi_file, gca_file, directory, url_file, ena_only):
    with SqliteDict(ncbi_file) as ncbi:
        base = Path(directory)
        for line in gca_file:
            row = json.loads(line)
            gca = determine_latest(ncbi, row['accession'])
            save_path = base / f"{row['upi']}.fa.gz"
            url = ftp_path(ncbi, gca)
            if url is None:
                if gca.startswith("GCF_"):
                    raise ValueError(f"Somehow missing a RefSeq genome: {row}")
                ena_info = dict(row)
                ena_info['accession'] = None
                json.dump(ena_info, ena_only)
                ena_only.write('\n')
            else:
                url_file.write(str(save_path))
                url_file.write('\n')
                url_file.write(url)
                url_file.write('\n')


if __name__ == '__main__':
    main()
