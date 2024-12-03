# -*- coding: utf-8 -*-

"""
Copyright [2009-2023] EMBL-European Bioinformatics Institute
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

import pytest

from rfamseq import wgs
from rfamseq.accession import Accession
from rfamseq.missing import Missing


def test_missing_knows_it_is_empty():
    assert bool(Missing.empty()) is False
    assert (not Missing) is False


def test_missing_knows_it_is_full():
    missing = Missing.empty()
    missing.add(Accession.build("NC_004448.1"))
    assert bool(missing) is True


def test_it_tracks_what_is_given():
    missing = Missing.empty()
    missing.add(wgs.WgsSequenceId.build("AAFI02000001"))
    missing.add(wgs.WgsPrefix.build("ABFI02"))
    missing.add(Accession.build("NC_004448.1"))
    assert missing == Missing(
        accessions={Accession.build("NC_004448.1")},
        wgs_sets={wgs.WgsPrefix.build("ABFI02")},
        wgs_sequences={wgs.WgsSequenceId.build("AAFI02000001")},
    )


def test_add_wgs_suppresses_members():
    missing = Missing.empty()
    missing.add(wgs.WgsSequenceId.build("AAFI02000001"))
    missing.add(wgs.WgsPrefix.build("AAFI02"))
    assert missing == Missing(
        accessions=set(), wgs_sequences=set(), wgs_sets={wgs.WgsPrefix.build("AAFI02")}
    )
