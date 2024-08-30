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

from rfamseq.ncbi import ftp


@pytest.mark.parametrize(
    "given,expected",
    [
        (
            "GCA_000292645.1",
            "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/292/645/GCA_000292645.1_ASM29264v1/GCA_000292645.1_ASM29264v1_genomic.fna.gz",
        ),
        (
            "GCA_000499845.1",
            "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/499/845/GCA_000499845.1_PhaVulg1_0/GCA_000499845.1_PhaVulg1_0_genomic.fna.gz",
        ),
        (
            "GCA_000499845",
            "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/499/845/GCA_000499845.1_PhaVulg1_0/GCA_000499845.1_PhaVulg1_0_genomic.fna.gz",
        ),
        (
            "GCA_000411955.5",
            "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/411/955/GCA_000411955.5_PG29_v4.1/GCA_000411955.5_PG29_v4.1_genomic.fna.gz",
        ),
        (
            "GCA_000208925.2",
            "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/208/925/GCA_000208925.2_JCVI_ESG2_1.0/GCA_000208925.2_JCVI_ESG2_1.0_genomic.fna.gz",
        ),
        ("GCA_000411955", None),
        ("000499845", None),
        ("", None),
    ],
)
def test_can_correctly_guess_ftp_url(given, expected):
    assert ftp.guess_genome_url(given, "_genomic.fna.gz") == expected
