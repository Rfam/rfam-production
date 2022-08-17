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

import pytest

from rfamseq import ncbi

# import tempfile
# from sqlitedict import SqliteDict


# @pytest.fixture(scope='module')
# def info() -> ty.Generator[SqliteDict]:
#     with tempfile.NamedTemporaryFile() as tmp:
#         with SqliteDict(tmp.name) as db:
#             yield db


@pytest.mark.parametrize(
    "path,expected",
    [
        (
            "tests/data/GCA_000411955.6_PG29_v5_assembly_report.txt",
            ncbi.NcbiAssemblyInfo(
                taxid=3330,
                assembly_name="PG29_v5",
                assembly_level=ncbi.AssemblyLevel.SCAFFOLD,
                organism_name="Picea glauca (white spruce)",
                bio_sample="SAMN01120252",
                bio_project="PRJNA83435",
                wgs_project="ALWZ05",
                sequence_info=[
                    ncbi.NcbiSequenceInfo(
                        genbank_accession="ALWZ05S0000001.1",
                        name="scaffoldPg-05r180724s0000001",
                        role=ncbi.SequenceRole.UNPLACED_SCAFFOLD,
                        molecule_type=None,
                        length=1717731,
                    ),
                    ncbi.NcbiSequenceInfo(
                        genbank_accession="ALWZ05S0000002.1",
                        name="scaffoldPg-05r180724s0000002",
                        role=ncbi.SequenceRole.UNPLACED_SCAFFOLD,
                        molecule_type=None,
                        length=2066984,
                    )
                ],
            ),
        ),
        (
            "tests/data/GCA_000597845.1_ASM59784v1_assembly_report.txt",
            ncbi.NcbiAssemblyInfo(
                taxid=562,
                assembly_name="ASM59784v1",
                assembly_level=ncbi.AssemblyLevel.COMPLETE_GENOME,
                organism_name="Escherichia coli (E. coli)",
                bio_sample="SAMN02647163",
                bio_project="PRJNA238952",
                wgs_project=None,
                sequence_info=[
                    ncbi.NcbiSequenceInfo(
                        genbank_accession="CP007265.1",
                        name="ANONYMOUS",
                        role=ncbi.SequenceRole.ASSEMBLED_MOLECULE,
                        molecule_type="Chromosome",
                        length=4758629,
                    )
                ],
            ),
        ),
        (
            "tests/data/GCA_001963235.1_ViralProj361950_assembly_report.txt",
            ncbi.NcbiAssemblyInfo(
                taxid=1919083,
                assembly_name="ViralProj361950",
                assembly_level=None,
                organism_name="Cynomolgus cytomegalovirus (viruses)",
                bio_sample=None,
                bio_project="PRJNA485481",
                wgs_project=None,
                sequence_info=[
                    ncbi.NcbiSequenceInfo(
                        genbank_accession="KX689265.1",
                        name="NC_033176.1",
                        role=ncbi.SequenceRole.ASSEMBLED_MOLECULE,
                        molecule_type="Segment",
                        length=None,
                    )
                ],
            ),
        ),
    ],
)
def test_can_parse_assemblies(path, expected):
    with open(path, "r") as raw:
        assert ncbi.parse_assembly_info(raw) == expected
