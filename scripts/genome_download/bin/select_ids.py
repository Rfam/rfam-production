#!/usr/bin/env python

import shutil
import json
from pathlib import Path
import typing as ty
import tempfile
import subprocess as sp
from io import StringIO

from Bio import SeqIO

import click

import download_ena


def versionless(record: SeqIO.SeqRecord) -> str:
    return record.id.split('.', 1)[0]


def select_ids(path: Path, allowed: ty.Set[str]) -> ty.Iterable[SeqIO.SeqRecord]:
    seen = set()
    with path.open('r') as raw:
        for record in SeqIO.parse(raw, 'fasta'):
            for key in [record.id, versionless(record)]:
                if key not in allowed:
                    continue
                if key in seen:
                    raise ValueError(f"Duplicate ids in fasta file at {path}")
                seen.add(key)
                yield record
                break

    if len(seen) != len(allowed):
        missing = allowed - seen
        yield from download_ena.fetch(missing)


def load_ignore(handle: ty.IO) -> ty.Set[str]:
    ignore = set()
    for line in handle:
        data = json.loads(line)
        ignore.add(data['upi'])
    return ignore


@click.command()
@click.option('--ignore-file', type=click.Path())
@click.argument('gca-file', type=click.File('r'))
@click.argument('directory', type=click.Path())
@click.argument('output', default='.', type=click.Path())
def main(gca_file: ty.IO, directory: str, output: str, ignore_file: ty.Optional[str]=None):
    base = Path(directory)
    out = Path(output)
    ignore = set()
    if ignore_file and Path(ignore_file).exists():
        with open(ignore_file, 'r') as raw:
            ignore = load_ignore(raw)
    for line in gca_file:
        row = json.loads(line)
        upi = row['upi']
        path = base / f'{upi}.fa'
        if not path.exists() and upi not in ignore:
            raise ValueError(f"Failed to find genome for {upi} at {path}")

        final_path = out / f'{upi}.fa'
        if row['ids'] == ['*']:
            shutil.copyfile(path, final_path)
            return

        allowed = {id for id in row['ids']}
        with open(final_path, 'w') as f:
            SeqIO.write(select_ids(path, allowed), f, 'fasta')


if __name__ == '__main__':
    main()
