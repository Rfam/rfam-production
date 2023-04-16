# -*- coding: utf-8 -*-

"""
Copyright [2009-${2022}] EMBL-European Bioinformatics Institute
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

import cattrs
import pytest

from rfamseq.accession import Accession


@pytest.mark.parametrize(
    "raw,expected",
    [
        ("NC_004448.1", Accession(accession="NC_004448", version="1")),
        ("NC_004448", Accession(accession="NC_004448", version=None)),
        ("CY130004", Accession(accession="CY130004", version=None)),
        ("D00050", Accession(accession="D00050", version=None)),
    ],
)
def test_can_build_as_expected(raw, expected):
    assert Accession.build(raw) == expected


@pytest.mark.parametrize(
    "raw,expected",
    [
        ("NC_004448.1", Accession(accession="NC_004448", version="1")),
        ("NC_004448", Accession(accession="NC_004448", version=None)),
        (
            {"accession": "NC_004448", "version": None},
            Accession(accession="NC_004448", version=None),
        ),
        (
            {"accession": "NC_004448", "version": "1"},
            Accession(accession="NC_004448", version="1"),
        ),
    ],
)
def test_can_be_structured_by_cattrs(raw, expected):
    assert cattrs.structure(raw, Accession) == expected


@pytest.mark.parametrize(
    "raw",
    [
        "NC_004448.1",
        "NC_004448",
        "CY130004",
        "D00050",
    ],
)
def test_can_round_trip_to_a_string(raw):
    assert str(Accession.build(raw)) == raw


@pytest.mark.parametrize(
    "raw",
    [
        "NC_004448.1",
        "NC_004448",
        "CY130004",
        "D00050",
    ],
)
def test_can_round_trip_with_cattrs_to_a_string(raw):
    assert str(cattrs.structure(raw, Accession)) == raw


@pytest.mark.parametrize(
    "given,other,expected",
    [
        (Accession.build("NC_004448.1"), Accession.build("NC_004448.1"), True),
        (Accession.build("NC_004448.1"), Accession.build("NC_004448.2"), True),
        (Accession.build("NC_004448.1"), Accession.build("NC_004448"), True),
        (Accession.build("NC_004448"), Accession.build("NC_004448.1"), True),
        (Accession.build("NC_004448"), Accession.build("NC_004448"), True),
        (Accession.build("NC_004449"), Accession.build("NC_004448"), False),
        (
            Accession.build("NC_004449.1", aliases=(Accession.build("NC_004448"),)),
            Accession.build("NC_004448.1"),
            True,
        ),
        (
            Accession.build("NC_004448.1"),
            Accession.build("NC_004448.1", aliases=(Accession.build("NC_004448"),)),
            True,
        ),
        (
            Accession.build("NC_004449", aliases=(Accession.build("NC_004448"),)),
            Accession.build("NC_004448.1"),
            True,
        ),
        (
            Accession.build("NC_004448"),
            Accession.build("NC_004448.1", aliases=(Accession.build("NC_004448"),)),
            True,
        ),
    ],
)
def test_can_detect_expected_matches(given, other, expected):
    assert given.matches(other) == expected
