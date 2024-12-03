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

import contextlib
import os
from pathlib import Path

import pytest

from rfamseq import ena, wgs


@contextlib.contextmanager
def set_env(**environ):
    """
    Temporarily set the process environment variables.

    >>> with set_env(PLUGINS_DIR='test/plugins'):
    ...   "PLUGINS_DIR" in os.environ
    True

    >>> "PLUGINS_DIR" in os.environ
    False

    :type environ: dict[str, unicode]
    :param environ: Environment variables to set
    """
    old_environ = dict(os.environ)
    os.environ.update(environ)
    try:
        yield
    finally:
        os.environ.clear()
        os.environ.update(old_environ)


@pytest.mark.parametrize(
    "accession,expected",
    [
        (
            "CABU01000000",
            "https://g-a8b222.dd271.03c0.data.globus.org/pub/databases/ena/wgs/public/cab/CABU01.fasta.gz",
        ),
        (
            "CABU01",
            "https://g-a8b222.dd271.03c0.data.globus.org/pub/databases/ena/wgs/public/cab/CABU01.fasta.gz",
        ),
        (
            "JABWAI01",
            "https://g-a8b222.dd271.03c0.data.globus.org/pub/databases/ena/wgs/public/jab/JABWAI01.fasta.gz",
        ),
    ],
)
def test_can_generate_expected_wgs_fasta_url(accession, expected):
    assert ena.wgs_fasta_url(wgs.WgsPrefix.build(accession)) == expected


@pytest.mark.parametrize(
    "accession,expected",
    [
        (
            "CABU01000000",
            "https://g-a8b222.dd271.03c0.data.globus.org/pub/databases/ena/wgs/suppressed/cab/CABU01.fasta.gz",
        ),
        (
            "CABU01",
            "https://g-a8b222.dd271.03c0.data.globus.org/pub/databases/ena/wgs/suppressed/cab/CABU01.fasta.gz",
        ),
        (
            "JABWAI01",
            "https://g-a8b222.dd271.03c0.data.globus.org/pub/databases/ena/wgs/suppressed/jab/JABWAI01.fasta.gz",
        ),
    ],
)
def test_can_generate_expected_suppressed_wgs_fasta_url(accession, expected):
    assert (
        ena.wgs_fasta_url(wgs.WgsPrefix.build(accession), use_suppressed=True)
        == expected
    )


@pytest.mark.parametrize(
    "given,expected",
    [
        (
            "ftp://ftp.ebi.ac.uk/pub/databases/ena/wgs/public/jab/JABWAI01.fasta.gz",
            Path("/ftp/ena/databases/ena/wgs/public/jab/JABWAI01.fasta.gz"),
        ),
        (
            "https://g-a8b222.dd271.03c0.data.globus.org/pub/databases/ena/wgs/public/jab/JABWAI01.fasta.gz",
            Path("/ftp/ena/databases/ena/wgs/public/jab/JABWAI01.fasta.gz"),
        ),
    ],
)
def test_can_convert_to_expected_internal_paths(given, expected):
    with set_env(ENA_PATH="/ftp/ena"):
        assert ena.internal_path(given) == expected
