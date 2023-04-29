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
import typing as ty
from contextlib import contextmanager
from pathlib import Path
from urllib.parse import urlparse

from rfamseq import wget, wgs

LOGGER = logging.getLogger(__name__)

ENA_FASTA_URL = "https://www.ebi.ac.uk/ena/browser/api/fasta/{accession}?download=true"
ENA_EMBL_URL = "https://www.ebi.ac.uk/ena/browser/api/embl/{accession}?download=true"
ENA_WGS_FASTA_URL = (
    "ftp://ftp.ebi.ac.uk/pub/databases/ena/wgs/public/{prefix}/{name}.fasta.gz"
)
ENA_SUPPRESED_WGS_FASTA_URL = (
    "ftp://ftp.ebi.ac.uk/pub/databases/ena/wgs/suppressed/{prefix}/{name}.fasta.gz"
)


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
        if not filepath.exists():
            LOGGER.info("Path %s does not exist", filepath)
        else:
            with filepath.open("r") as handle:
                LOGGER.debug("Using local file path %s", filepath)
                yield handle
            return

    with wget.wget(url) as handle:
        LOGGER.debug("Using FTP fetch of %s", url)
        yield handle


@contextmanager
def fetch_fasta(accession: str) -> ty.Iterator[ty.IO]:
    LOGGER.info("Fetching %s fasta from ENA", accession)
    with fetch(ENA_FASTA_URL, accession=accession) as handle:
        yield handle


def wgs_fasta_url(prefix: wgs.WgsPrefix, use_suppressed=False) -> str:
    short = prefix.wgs_id[0:3].lower()
    name = prefix.to_wgs_string().upper()
    name = re.sub("0+$", "", name)
    if use_suppressed:
        return ENA_SUPPRESED_WGS_FASTA_URL.format(prefix=short, name=name)
    return ENA_WGS_FASTA_URL.format(prefix=short, name=name)


@contextmanager
def wgs_fasta(
    prefix: wgs.WgsPrefix, max_increase=2, use_suppressed=False
) -> ty.Iterator[ty.IO]:
    url = wgs_fasta_url(prefix, use_suppressed=use_suppressed)
    LOGGER.info("Fetching the wgs fasta set for %s at %s", prefix, url)
    try:
        with fetch(url) as handle:
            yield handle
    except wget.FetchError as err:
        LOGGER.exception(err)
        if max_increase <= 0:
            raise err
        LOGGER.info("Trying to increment and grab next WGS set")
        next_url = wgs_fasta_url(prefix.next_version())
        with fetch(next_url, max_increase=max_increase - 1) as handle:
            yield handle
