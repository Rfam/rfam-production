# -*- coding: utf-8 -*-

"""
Copyright [2009-2022] EMBL-European Bioinformatics Institute
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

import logging
import os
import re
import tempfile
import typing as ty
from contextlib import contextmanager
from pathlib import Path
from urllib.parse import urlparse

from Bio import SeqIO

from rfamseq import fasta, wget, wgs, rsync

LOGGER = logging.getLogger(__name__)

ENA_FASTA_URL = "https://www.ebi.ac.uk/ena/browser/api/fasta/{accession}?download=true"
ENA_EMBL_URL = "https://www.ebi.ac.uk/ena/browser/api/embl/{accession}?download=true"
ENA_WGS_FASTA_URL = (
    "ftp://ftp.ebi.ac.uk/pub/databases/ena/wgs/public/{prefix}/{name}.fasta.gz"
)
ENA_SUPPRESED_WGS_FASTA_URL = (
    "ftp://ftp.ebi.ac.uk/pub/databases/ena/wgs/suppressed/{prefix}/{name}.fasta.gz"
)

GLOBUS_WGS_FASTA_URL = "https://g-a8b222.dd271.03c0.data.globus.org/pub/databases/ena/wgs/public/{prefix}/{name}.fasta.gz"
GLOBUS_SUPPRESED_WGS_FASTA_URL = "https://g-a8b222.dd271.03c0.data.globus.org/pub/databases/ena/wgs/suppressed/{prefix}/{name}.fasta.gz"


def internal_path(url: str) -> ty.Optional[Path]:
    base_path = os.environ.get("ENA_PATH", None)
    if not base_path:
        return None
    parsed = urlparse(url)
    if "ftp" not in parsed.netloc:
        return None
    parts = "/".join(parsed.path.split("/")[2:])
    path = Path(base_path)
    path = path / Path(parts)
    return path


@contextmanager
def fetch(template: str, **data) -> ty.Iterator[ty.IO]:
    url = template.format(**data)
    LOGGER.debug("Fetching %s", url)
    if filepath := internal_path(url):
        LOGGER.info("Trying to use internal path %s", filepath)
        try:
            with rsync.rsync(filepath) as handle:
                LOGGER.debug("Using local file path %s", filepath)
                yield handle
            return
        except Exception:
            LOGGER.info("Could not open %s", filepath)

    with wget.wget(url) as handle:
        LOGGER.debug("Using FTP fetch of %s", url)
        yield handle


@contextmanager
def normalized(template: str, **data) -> ty.Iterator[ty.IO]:
    with tempfile.NamedTemporaryFile(mode="w", dir=os.curdir) as tmp:
        with fetch(template, **data) as raw:
            SeqIO.write(fasta.parse(raw), tmp.name, "fasta")
        tmp.flush()
        tmp.seek(0)
        yield tmp


@contextmanager
def fetch_fasta(accession: str) -> ty.Iterator[ty.IO]:
    LOGGER.info("Fetching %s fasta from ENA", accession)
    with normalized(ENA_FASTA_URL, accession=accession) as handle:
        yield handle


def wgs_fasta_url(prefix: wgs.WgsPrefix, use_suppressed=False) -> str:
    short = prefix.wgs_id[0:3].lower()
    name = prefix.to_wgs_string().upper()
    name = re.sub("0+$", "", name)
    if use_suppressed:
        return GLOBUS_SUPPRESED_WGS_FASTA_URL.format(prefix=short, name=name)
    return GLOBUS_WGS_FASTA_URL.format(prefix=short, name=name)


@contextmanager
def wgs_fasta(
    prefix: wgs.WgsPrefix, max_increase=2, use_suppressed=False
) -> ty.Iterator[ty.IO]:
    url = wgs_fasta_url(prefix, use_suppressed=use_suppressed)
    LOGGER.info("Fetching the wgs fasta set for %s at %s", prefix, url)
    try:
        with normalized(url) as handle:
            yield handle
    except Exception as err:
        LOGGER.exception(err)
        if max_increase <= 0:
            raise err
        LOGGER.info("Trying to increment and grab next WGS set")
        next_url = wgs_fasta_url(prefix.next_version())
        with normalized(next_url, max_increase=max_increase - 1) as handle:
            yield handle
