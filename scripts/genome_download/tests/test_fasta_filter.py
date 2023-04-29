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

from pathlib import Path

import pytest

from rfamseq import ncbi, uniprot
from rfamseq.accession import Accession
from rfamseq.fasta_filter import FastaFilter
from rfamseq.wgs import WgsPrefix, WgsSequenceId


@pytest.mark.parametrize(
    "upi,genome,given,expected",
    [
        (
            "UP000007062",
            "GCA_000005575.1_AgamP3_assembly_report.txt",
            Accession.build("CM000360.1"),
            [Accession.build("CM000360")],
        ),
        (
            "UP000007062",
            "GCA_000005575.1_AgamP3_assembly_report.txt",
            Accession.build("CM000360"),
            [Accession.build("CM000360")],
        ),
        (
            "UP000007062",
            "GCA_000005575.1_AgamP3_assembly_report.txt",
            Accession.build("NT_078266.2"),
            [Accession.build("CM000357")],
        ),
        (
            "UP000007062",
            "GCA_000005575.1_AgamP3_assembly_report.txt",
            Accession.build("NT_078266"),
            [Accession.build("CM000357")],
        ),
        (
            "UP000007062",
            "GCA_000005575.1_AgamP3_assembly_report.txt",
            Accession.build("BM000999"),
            [],
        ),
        (
            "UP000007062",
            "GCA_000005575.1_AgamP3_assembly_report.txt",
            WgsSequenceId.build("AAAB01001191.1"),
            [WgsPrefix.build("AAAB01")],
        ),
        # (
        #     "UP000007062",
        #     "GCA_000005575.1_AgamP3_assembly_report.txt",
        #     Accession.build("NW_164059.1"),
        #     [WgsPrefix.build("AAAB01")],
        # ),
    ],
)
def test_can_detect_expected_matches(upi, genome, given, expected):
    assembly_report = None
    genome_path = Path(f"tests/data/{genome}")
    if genome_path.exists():
        with genome_path.open("r") as raw:
            assembly_report = ncbi.parse_assembly_info(raw)
            assert assembly_report

    proteome_path = Path(f"tests/data/{upi}.xml")
    info = list(uniprot.proteomes(proteome_path, set()))
    assert len(info) == 1
    info = info[0]
    requested = info.genome_info.components
    filter = FastaFilter.from_selected(assembly_report, requested)
    assert filter.matching_components(given) == expected
