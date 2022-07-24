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

import typing as ty

import tempfile
from sqlitedict import SqliteDict

import pytest

from rfamseq import ncbi


@pytest.fixture(scope='module')
def info() -> ty.Generator[SqliteDict]:
    with tempfile.NamedTemporaryFile() as tmp:
        with SqliteDict(tmp.name) as db:
            yield db


@pytest.paramterize('accession,expected', [
    ()
])
def test_ftp_path(info, accession, expected):
    assert ncbi.ftp_path(info, accession) == expected
