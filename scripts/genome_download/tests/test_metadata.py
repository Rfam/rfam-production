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

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from rfamseq import metadata, ncbi, uniprot


@pytest.mark.parametrize(
    "given,expected",
    [
        (ncbi.AssemblyLevel.CHROMOSOME, metadata.AssemblyLevel.CHROMOSOME),
        (ncbi.AssemblyLevel.SCAFFOLD, metadata.AssemblyLevel.SCAFFOLD),
        (ncbi.AssemblyLevel.CONTIG, metadata.AssemblyLevel.CONTIG),
    ],
)
def test_can_convert_assembly_levels(given, expected):
    assert metadata.AssemblyLevel.from_ncbi(given) == expected


@pytest.mark.parametrize(
    "given,expected",
    [
        (
            uniprot.LineageInfo(
                ncbi_id=10090,
                species="Mus musculus",
                common_name="Mouse",
                tax_string="Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Glires; Rodentia; Myomorpha; Muroidea; Muridae; Murinae; Mus; Mus.",
            ),
            metadata.Taxonomy(
                ncbi_id=10090,
                species="Mus musculus (mouse)",
                tax_string="Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Glires; Rodentia; Myomorpha; Muroidea; Muridae; Murinae; Mus; Mus.",
                tree_display_name="Mus_musculus_(mouse)",
                align_display_name="Mus_musculus_(mouse)[10090]",
            ),
        ),
        (
            uniprot.LineageInfo(
                ncbi_id=562,
                species="Escherichia coli",
                common_name=None,
                tax_string="Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacterales; Enterobacteriaceae; Escherichia.",
            ),
            metadata.Taxonomy(
                ncbi_id=562,
                species="Escherichia coli",
                tax_string="Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacterales; Enterobacteriaceae; Escherichia.",
                tree_display_name="Escherichia_coli",
                align_display_name="Escherichia_coli[562]",
            ),
        ),
        (
            uniprot.LineageInfo(
                ncbi_id=652674,
                species="Hepatitis E virus genotype 1 (isolate Human/China/HeBei/1987)",
                common_name="HEV",
                tax_string="Viruses; Riboviria; Orthornavirae; Kitrinoviricota; Alsuviricetes; Hepelivirales; Hepeviridae; Orthohepevirus; Hepatitis E virus.",
            ),
            metadata.Taxonomy(
                ncbi_id=652674,
                species="Hepatitis E virus genotype 1 (isolate Human/China/HeBei/1987)",
                tax_string="Viruses; Riboviria; Orthornavirae; Kitrinoviricota; Alsuviricetes; Hepelivirales; Hepeviridae; Orthohepevirus; Hepatitis E virus.",
                tree_display_name="Hepatitis_E_virus_genotype_1_(isolate_Human/China/HeBei/1987)",
                align_display_name="Hepatitis_E_virus_genotype_1_(isolate_Human/China/HeBei/1987)[652674]",
            ),
        ),
    ],
)
def test_can_build_taxonomy_from_lineage(given, expected):
    assert metadata.Taxonomy.from_lineage(given) == expected


@pytest.mark.parametrize(
    "given,expected",
    [
        (
            SeqRecord(Seq("aaa"), id="something", description="Hello world"),
            metadata.FromFasta(
                rfamseq_acc="something", length=len("aaa"), description="Hello world"
            ),
        ),
    ],
)
def test_can_build_from_fasta_entries(given, expected):
    assert metadata.FromFasta.from_record(given) == expected


@pytest.mark.parametrize(
    "taxid,info,expected",
    [
        (
            10,
            metadata.FromFasta(rfamseq_acc="ABC", length=10, description="sequence"),
            metadata.RfamSeq(
                rfamseq_acc="ABC",
                accession="ABC",
                version="000001",
                ncbi_id=10,
                mol_type=metadata.MoleculeType.GENOMIC_DNA,
                length=10,
                description="sequence",
                previous_acc="",
                source="UNIPROT; ENA",
            ),
        ),
        (
            10,
            metadata.FromFasta(
                rfamseq_acc="AAAA02006935.1", length=10, description="sequence"
            ),
            metadata.RfamSeq(
                rfamseq_acc="AAAA02006935.1",
                accession="AAAA02006935",
                version="000001",
                ncbi_id=10,
                mol_type=metadata.MoleculeType.GENOMIC_DNA,
                length=10,
                description="sequence",
                previous_acc="",
                source="UNIPROT; ENA",
            ),
        ),
        (
            10,
            metadata.FromFasta(
                rfamseq_acc="AACE03000008.2", length=10, description="sequence"
            ),
            metadata.RfamSeq(
                rfamseq_acc="AACE03000008.2",
                accession="AACE03000008",
                version="000002",
                ncbi_id=10,
                mol_type=metadata.MoleculeType.GENOMIC_DNA,
                length=10,
                description="sequence",
                previous_acc="",
                source="UNIPROT; ENA",
            ),
        ),
    ],
)
def test_can_build_rfamseq_entries(taxid, info, expected):
    assert metadata.RfamSeq.from_fasta(taxid, info) == expected


@pytest.mark.parametrize(
    "upid,version,info,expected",
    [
        (
            "UP1",
            "15.0",
            metadata.FromFasta(
                rfamseq_acc="AACE03000008.2", length=10, description="sequence"
            ),
            metadata.GenSeq(
                upid="UP1",
                rfamseq_acc="AACE03000008.2",
                chromosome_name=None,
                chromosome_type=None,
                version="15.0",
            ),
        ),
    ],
)
def test_can_build_genseq_from_entries(upid, version, info, expected):
    assert metadata.GenSeq.from_fasta(upid, version, info) == expected
