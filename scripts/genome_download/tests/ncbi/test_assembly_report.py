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
from io import StringIO

import pytest
import requests

from rfamseq import ncbi
from rfamseq.accession import Accession
from rfamseq.ncbi import assembly_report as report


@pytest.mark.parametrize(
    "path,expected",
    [
        (
            "tests/data/GCA_000411955.6_PG29_v5_assembly_report.txt",
            ncbi.NcbiAssemblyReport(
                taxid=3330,
                assembly_name="PG29_v5",
                assembly_level=ncbi.AssemblyLevel.SCAFFOLD,
                organism_name="Picea glauca (white spruce)",
                bio_sample="SAMN01120252",
                bio_project="PRJNA83435",
                wgs_project="ALWZ05",
                sequence_info=[
                    ncbi.NcbiSequenceInfo(
                        genbank_accession=Accession.build("ALWZ05S0000001.1"),
                        refseq_accession=None,
                        relationship=ncbi.SequenceRelationship.DIFFERENT,
                        name="scaffoldPg-05r180724s0000001",
                        role=ncbi.SequenceRole.UNPLACED_SCAFFOLD,
                        molecule_type=None,
                        length=1717731,
                    ),
                    ncbi.NcbiSequenceInfo(
                        genbank_accession=Accession.build("ALWZ05S0000002.1"),
                        refseq_accession=None,
                        relationship=ncbi.SequenceRelationship.DIFFERENT,
                        name="scaffoldPg-05r180724s0000002",
                        role=ncbi.SequenceRole.UNPLACED_SCAFFOLD,
                        molecule_type=None,
                        length=2066984,
                    ),
                ],
            ),
        ),
        (
            "tests/data/GCA_000597845.1_ASM59784v1_assembly_report.txt",
            ncbi.NcbiAssemblyReport(
                taxid=562,
                assembly_name="ASM59784v1",
                assembly_level=ncbi.AssemblyLevel.COMPLETE_GENOME,
                organism_name="Escherichia coli (E. coli)",
                bio_sample="SAMN02647163",
                bio_project="PRJNA238952",
                wgs_project=None,
                sequence_info=[
                    ncbi.NcbiSequenceInfo(
                        genbank_accession=Accession.build("CP007265.1"),
                        refseq_accession=Accession.build("NZ_CP007265.1"),
                        relationship=ncbi.SequenceRelationship.EQUAL,
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
            ncbi.NcbiAssemblyReport(
                taxid=1919083,
                assembly_name="ViralProj361950",
                assembly_level=None,
                organism_name="Cynomolgus cytomegalovirus (viruses)",
                bio_sample=None,
                bio_project="PRJNA485481",
                wgs_project=None,
                sequence_info=[
                    ncbi.NcbiSequenceInfo(
                        genbank_accession=Accession.build("KX689265.1"),
                        refseq_accession=Accession.build("NC_033176.1"),
                        relationship=ncbi.SequenceRelationship.EQUAL,
                        name="NC_033176.1",
                        role=ncbi.SequenceRole.ASSEMBLED_MOLECULE,
                        molecule_type="Segment",
                        length=None,
                    )
                ],
            ),
        ),
        (
            "tests/data/GCA_900105925.1_IMG-taxon_2639762615_annotated_assembly_assembly_report.txt",
            ncbi.NcbiAssemblyReport(
                taxid=419479,
                assembly_name="IMG-taxon 2639762615 annotated assembly",
                assembly_level=ncbi.AssemblyLevel.CHROMOSOME,
                organism_name="Jiangella alkaliphila (high GC Gram+)",
                bio_sample="SAMN04488563",
                bio_project="PRJEB16466",
                wgs_project=None,
                sequence_info=[
                    ncbi.NcbiSequenceInfo(
                        genbank_accession=Accession.build("LT629791.1"),
                        refseq_accession=Accession.build("NZ_LT629791.1"),
                        relationship=ncbi.SequenceRelationship.EQUAL,
                        name="I",
                        role=ncbi.SequenceRole.ASSEMBLED_MOLECULE,
                        molecule_type="Chromosome",
                        length=7716600,
                    )
                ],
            ),
        ),
        (
            "tests/data/GCF_000455745.1_ASM45574v1_assembly_report.txt",
            ncbi.NcbiAssemblyReport(
                taxid=38654,
                assembly_name="ASM45574v1",
                assembly_level=ncbi.AssemblyLevel.SCAFFOLD,
                organism_name="Alligator sinensis (Chinese alligator)",
                bio_sample="SAMN02981559",
                bio_project="PRJNA215016",
                wgs_project="AVPB01",
                sequence_info=[
                    ncbi.NcbiSequenceInfo(
                        genbank_accession=Accession.build("KE698715.1"),
                        refseq_accession=Accession.build("NW_005847012.1"),
                        relationship=ncbi.SequenceRelationship.EQUAL,
                        name="scaffold3284_1",
                        role=ncbi.SequenceRole.UNPLACED_SCAFFOLD,
                        molecule_type=None,
                        length=500,
                    ),
                    ncbi.NcbiSequenceInfo(
                        genbank_accession=None,
                        refseq_accession=Accession.build("NC_004448.1"),
                        relationship=ncbi.SequenceRelationship.DIFFERENT,
                        name="MT",
                        role=ncbi.SequenceRole.ASSEMBLED_MOLECULE,
                        molecule_type="Mitochondrion",
                        length=16746,
                    ),
                ],
            ),
        ),
    ],
)
def test_can_parse_assemblies(path, expected):
    with open(path, "r") as raw:
        assert ncbi.parse_assembly_info(raw) == expected


@pytest.mark.parametrize(
    "info,expected",
    [
        (
            ncbi.NcbiSequenceInfo(
                genbank_accession=None,
                refseq_accession=Accession.build("NC_004448.1"),
                relationship=ncbi.SequenceRelationship.DIFFERENT,
                name="MT",
                role=ncbi.SequenceRole.ASSEMBLED_MOLECULE,
                molecule_type="Mitochondrion",
                length=16746,
            ),
            Accession.build("NC_004448.1"),
        ),
        (
            ncbi.NcbiSequenceInfo(
                genbank_accession=Accession.build("KE698715.1"),
                refseq_accession=Accession.build("NW_005847012.1"),
                relationship=ncbi.SequenceRelationship.EQUAL,
                name="scaffold3284_1",
                role=ncbi.SequenceRole.UNPLACED_SCAFFOLD,
                molecule_type=None,
                length=500,
            ),
            Accession.build("KE698715.1", aliases=(Accession.build("NW_005847012.1"),)),
        ),
        (
            ncbi.NcbiSequenceInfo(
                genbank_accession=Accession.build("KE698715.1"),
                refseq_accession=Accession.build("NW_005847012.1"),
                relationship=ncbi.SequenceRelationship.DIFFERENT,
                name="scaffold3284_1",
                role=ncbi.SequenceRole.UNPLACED_SCAFFOLD,
                molecule_type=None,
                length=500,
            ),
            Accession.build("KE698715.1"),
        ),
        (
            ncbi.NcbiSequenceInfo(
                genbank_accession=Accession.build("ALWZ05S0000002.1"),
                refseq_accession=None,
                relationship=ncbi.SequenceRelationship.DIFFERENT,
                name="scaffoldPg-05r180724s0000002",
                role=ncbi.SequenceRole.UNPLACED_SCAFFOLD,
                molecule_type=None,
                length=2066984,
            ),
            Accession.build("ALWZ05S0000002.1"),
        ),
    ],
)
def test_can_build_expected_accessions(info, expected):
    assert info.accession() == expected


@pytest.mark.parametrize(
    "info,expected",
    [
        (
            ncbi.NcbiSequenceInfo(
                genbank_accession=None,
                refseq_accession=Accession.build("NC_004448.1"),
                relationship=ncbi.SequenceRelationship.DIFFERENT,
                name="MT",
                role=ncbi.SequenceRole.ASSEMBLED_MOLECULE,
                molecule_type="Mitochondrion",
                length=16746,
            ),
            Accession.build("NC_004448.1"),
        ),
        (
            ncbi.NcbiSequenceInfo(
                genbank_accession=Accession.build("KE698715.1"),
                refseq_accession=Accession.build("NW_005847012.1"),
                relationship=ncbi.SequenceRelationship.EQUAL,
                name="scaffold3284_1",
                role=ncbi.SequenceRole.UNPLACED_SCAFFOLD,
                molecule_type=None,
                length=500,
            ),
            Accession.build("KE698715.1", aliases=(Accession.build("NW_005847012.1"),)),
        ),
        (
            ncbi.NcbiSequenceInfo(
                genbank_accession=Accession.build("KE698715.1"),
                refseq_accession=Accession.build("NW_005847012.1"),
                relationship=ncbi.SequenceRelationship.DIFFERENT,
                name="scaffold3284_1",
                role=ncbi.SequenceRole.UNPLACED_SCAFFOLD,
                molecule_type=None,
                length=500,
            ),
            Accession.build("KE698715.1", aliases=(Accession.build("NW_005847012.1"),)),
        ),
        (
            ncbi.NcbiSequenceInfo(
                genbank_accession=Accession.build("ALWZ05S0000002.1"),
                refseq_accession=None,
                relationship=ncbi.SequenceRelationship.DIFFERENT,
                name="scaffoldPg-05r180724s0000002",
                role=ncbi.SequenceRole.UNPLACED_SCAFFOLD,
                molecule_type=None,
                length=2066984,
            ),
            Accession.build("ALWZ05S0000002.1"),
        ),
    ],
)
def test_can_build_expected_merged_accessions(info, expected):
    assert info.accession(merged=True) == expected


@pytest.mark.parametrize(
    "info,accession,expected",
    [
        (
            ncbi.NcbiSequenceInfo(
                genbank_accession=Accession(accession="KE698606", version="1"),
                refseq_accession=Accession(accession="NW_005844920", version="1"),
                relationship=ncbi.SequenceRelationship.EQUAL,
                name="scaffold2549_1",
                role=ncbi.SequenceRole.UNPLACED_SCAFFOLD,
                molecule_type=None,
                length=2014,
            ),
            Accession.build("NC_004448.1"),
            False,
        ),
        (
            ncbi.NcbiSequenceInfo(
                genbank_accession=Accession(accession="KE698606", version="1"),
                refseq_accession=Accession(accession="NW_005844920", version="1"),
                relationship=ncbi.SequenceRelationship.EQUAL,
                name="scaffold2549_1",
                role=ncbi.SequenceRole.UNPLACED_SCAFFOLD,
                molecule_type=None,
                length=2014,
            ),
            Accession.build("NW_005844920"),
            True,
        ),
        (
            ncbi.NcbiSequenceInfo(
                genbank_accession=Accession(accession="KE698606", version="1"),
                refseq_accession=Accession(accession="NW_005844920", version="1"),
                relationship=ncbi.SequenceRelationship.EQUAL,
                name="scaffold2549_1",
                role=ncbi.SequenceRole.UNPLACED_SCAFFOLD,
                molecule_type=None,
                length=2014,
            ),
            Accession.build("NW_005844920.1"),
            True,
        ),
        (
            ncbi.NcbiSequenceInfo(
                genbank_accession=Accession(accession="KE698606", version="1"),
                refseq_accession=Accession(accession="NW_005844920", version="1"),
                relationship=ncbi.SequenceRelationship.EQUAL,
                name="scaffold2549_1",
                role=ncbi.SequenceRole.UNPLACED_SCAFFOLD,
                molecule_type=None,
                length=2014,
            ),
            Accession.build("KE698606.1"),
            True,
        ),
        (
            ncbi.NcbiSequenceInfo(
                genbank_accession=Accession(accession="KE698606", version="1"),
                refseq_accession=Accession(accession="NW_005844920", version="1"),
                relationship=ncbi.SequenceRelationship.EQUAL,
                name="scaffold2549_1",
                role=ncbi.SequenceRole.UNPLACED_SCAFFOLD,
                molecule_type=None,
                length=2014,
            ),
            Accession.build("KE698606"),
            True,
        ),
    ],
)
def test_sequence_info_detects_matches(info, accession, expected):
    assert info.matches(accession) == expected


@pytest.mark.parametrize(
    "url,accession,expected",
    [
        (
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/455/745/GCF_000455745.1_ASM45574v1/GCF_000455745.1_ASM45574v1_assembly_report.txt",
            "NW_005841829.1",
            True,
        )
    ],
)
def test_can_detect_when_sequence_is_unplaced(url, accession, expected):
    response = requests.get(url)
    response.raise_for_status()
    rep = report.parse_assembly_info(StringIO(response.text))
    assert rep
    assert rep.is_unplaced(Accession.build(accession)) == expected
