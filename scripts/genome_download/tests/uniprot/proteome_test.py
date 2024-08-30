# -*- coding: utf-8 -*-

"""
Copyright [2009-${2024}] EMBL-European Bioinformatics Institute
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

import json
import tempfile
from pathlib import Path

import pytest
import requests

from rfamseq.uniprot import proteome


def test_can_parse_current_reference():
    with tempfile.NamedTemporaryFile("w") as tmp:
        response = requests.get(
            "https://rest.uniprot.org/proteomes/stream?compressed=false&format=json&query=%28*%29+AND+%28proteome_type%3A1%29"
        )
        json.dump(response.json(), tmp)
        tmp.flush()

        data = list(proteome.parse(Path(tmp.name)))
    assert data


@pytest.mark.parametrize(
    "filename,count",
    [
        ("tests/data/proteomes.json", 24074),
    ],
)
def test_can_parse_expected_counts(filename: str, count: int):
    assert len(list(proteome.parse(Path(filename)))) == count


@pytest.mark.parametrize(
    "filename,expected",
    [
        ("tests/data/UP000214656.json", -1),
        ("tests/data/UP000192746.json", -1),
        ("tests/data/UP000177710.json", -1),
        ("tests/data/UP000054324.json", -1),
        ("tests/data/UP000054308.json", -1),
        ("tests/data/UP000214374.json", -1),
        ("tests/data/UP000189701.json", -1),
        ("tests/data/UP000178744.json", -1),
        ("tests/data/UP000186774.json", -1),
        ("tests/data/UP000007078.json", -1),
        ("tests/data/UP000094542.json", -1),
        ("tests/data/UP000007062.json", -1),
        ("tests/data/UP000214715.json", -1),
        ("tests/data/UP000007648.json", -1),
        ("tests/data/UP000000230.json", -1),
        ("tests/data/UP000186312.json", -1),
        ("tests/data/UP000000226.json", -1),
        ("tests/data/UP000214542.json", -1),
        ("tests/data/UP000006672.json", -1),
        ("tests/data/UP000214372.json", -1),
        ("tests/data/UP000008177.json", -1),
        ("tests/data/UP000011116.json", -1),
        ("tests/data/UP000000212.json", -1),
        ("tests/data/UP000093657.json", -1),
        ("tests/data/UP000229442.json", -1),
        ("tests/data/UP000000224.json", -1),
        ("tests/data/UP000001819.json", -1),
        ("tests/data/UP000033569.json", -1),
        ("tests/data/UP000006719.json", -1),
        ("tests/data/UP000019116.json", -1),
        ("tests/data/UP000005640.json", -1),
        ("tests/data/UP000198211.json", -1),
        ("tests/data/UP000000408.json", -1),
        ("tests/data/UP000012065.json", -1),
        ("tests/data/UP000662757.json", -1),
        ("tests/data/UP000177381.json", -1),
        ("tests/data/UP000035681.json", -1),
        ("tests/data/UP000214863.json", -1),
        ("tests/data/UP000214376.json", -1),
        ("tests/data/UP000006705.json", -1),
        ("tests/data/UP000214710.json", -1),
        ("tests/data/UP000035680.json", -1),
        ("tests/data/UP000189670.json", -1),
    ],
)
def test_can_parse_as_expected(filename: str, expected: proteome.Proteome):
    data = list(proteome.parse(Path(filename)))
    assert len(data) == 1
    assert data == expected


@pytest.mark.parametrize("expected", [])
def test_can_fetch_data(expected: proteome.Proteome):
    assert proteome.fetch(set([expected.id])) == expected
