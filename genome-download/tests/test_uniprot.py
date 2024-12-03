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
from pathlib import Path

import pytest

from rfamseq import uniprot
from rfamseq.accession import Accession
from rfamseq.wgs import WgsPrefix, WgsSequenceId

# UP000011116
# UP000012065
# UP000033569
# UP000054308
# UP000054324
# UP000094542
# UP000177381
# UP000177710
# UP000178744
# UP000186312
# UP000186774
# UP000189670
# UP000189701
# UP000192746
# UP000198211
# UP000214372
# UP000214374
# UP000214376


@pytest.mark.parametrize(
    "taxid,expected",
    [
        (
            "562",
            uniprot.LineageInfo(
                ncbi_id=562,
                species="Escherichia coli",
                common_name=None,
                tax_string="Bacteria; Pseudomonadota; Gammaproteobacteria; Enterobacterales; Enterobacteriaceae; Escherichia.",
            ),
        ),
        (
            "9606",
            uniprot.LineageInfo(
                ncbi_id=9606,
                species="Homo sapiens",
                common_name="Human",
                tax_string="Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini; Catarrhini; Hominidae; Homo.",
            ),
        ),
        (
            "10090",
            uniprot.LineageInfo(
                ncbi_id=10090,
                species="Mus musculus",
                common_name="Mouse",
                tax_string="Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Glires; Rodentia; Myomorpha; Muroidea; Muridae; Murinae; Mus; Mus.",
            ),
        ),
    ],
)
def test_can_get_correct_lineage(taxid, expected):
    assert uniprot.lineage_info(taxid) == expected


@pytest.mark.parametrize(
    "path,expected",
    [
        (
            Path("tests/data/UP000005640.xml"),
            uniprot.ProteomeInfo(
                upi="UP000005640",
                taxid=9606,
                is_reference=True,
                is_representative=True,
                proteome_description="""Homo sapiens (Homo sapiens sapiens) or modern humans are the only living species of the evolutionary branch of great apes known as hominids. Divergence of early humans from chimpanzees and gorillas is estimated to have occurred between 4 and 8 million years ago. The genus Homo (Homo habilis) appeared in Africa around 2.3 million years ago and shows the first signs of stone tool usage. The exact lineage of Homo species ie: H. habilis/H. ergaster to H. erectus to H. rhodesiensis/H.heidelbergensis to H. sapiens is still hotly disputed. However, continuing evolution and in particular larger brain size and complexity culminates in Homo sapiens. The first anatomically modern humans appear in the fossil record around 200,000 years ago. Modern humans migrated across the globe essentially as hunter-gatherers until around 12,000 years ago when the practice of agriculture and animal domestication enabled large populations to grow leading to the development of civilizations.""",
                genome_info=uniprot.GenomeInfo(
                    accession="GCA_000001405.27",
                    description=None,
                    source=uniprot.GenomeSource.ENSEMBL,
                    components=uniprot.SelectedComponents.build(
                        [
                            Accession.build("CM000663"),
                            Accession.build("CM000664"),
                            Accession.build("CM000665"),
                            Accession.build("CM000666"),
                            Accession.build("CM000667"),
                            Accession.build("CM000668"),
                            Accession.build("CM000669"),
                            Accession.build("CM000670"),
                            Accession.build("CM000671"),
                            Accession.build("CM000672"),
                            Accession.build("CM000673"),
                            Accession.build("CM000674"),
                            Accession.build("CM000675"),
                            Accession.build("CM000676"),
                            Accession.build("CM000677"),
                            Accession.build("CM000678"),
                            Accession.build("CM000679"),
                            Accession.build("CM000680"),
                            Accession.build("CM000681"),
                            Accession.build("CM000682"),
                            Accession.build("CM000683"),
                            Accession.build("CM000684"),
                            Accession.build("CM000685"),
                            Accession.build("CM000686"),
                            Accession.build("J01415"),
                            uniprot.UNPLACED,
                        ]
                    ),
                ),
                lineage_info=uniprot.LineageInfo(
                    ncbi_id=9606,
                    species="Homo sapiens",
                    common_name="Human",
                    tax_string="Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini; Catarrhini; Hominidae; Homo.",
                ),
            ),
        ),
        (
            Path("tests/data/UP000000226.xml"),
            uniprot.ProteomeInfo(
                upi="UP000000226",
                taxid=3885,
                is_reference=False,
                is_representative=True,
                proteome_description="The common bean is most widely cultivated of all beans in temperate regions, and is widely cultivated in semitropical regions. In some regions of the world dry beans furnish a large portion of the protein needs of low and middle class families.",
                genome_info=uniprot.GenomeInfo(
                    accession="GCA_000499845.1",
                    description=None,
                    source=uniprot.GenomeSource.ENA,
                    components=uniprot.SelectedComponents.build(
                        [
                            Accession.build("CM002288"),
                            Accession.build("CM002289"),
                            Accession.build("CM002290"),
                            Accession.build("CM002291"),
                            Accession.build("CM002292"),
                            Accession.build("CM002293"),
                            Accession.build("CM002294"),
                            Accession.build("CM002295"),
                            Accession.build("CM002296"),
                            Accession.build("CM002297"),
                            Accession.build("CM002298"),
                            WgsPrefix.build("ANNZ01000000"),
                        ]
                    ),
                ),
                lineage_info=uniprot.LineageInfo(
                    ncbi_id=3885,
                    species="Phaseolus vulgaris",
                    common_name="Kidney bean",
                    tax_string="Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta; Spermatophyta; Magnoliopsida; eudicotyledons; Gunneridae; Pentapetalae; rosids; fabids; Fabales; Fabaceae; Papilionoideae; 50 kb inversion clade; NPAAA clade; indigoferoid/millettioid clade; Phaseoleae; Phaseolus.",
                ),
            ),
        ),
        (
            Path("tests/data/UP000006672.xml"),
            uniprot.ProteomeInfo(
                upi="UP000006672",
                taxid=6279,
                is_reference=True,
                is_representative=True,
                proteome_description=None,
                genome_info=uniprot.GenomeInfo(
                    accession="GCA_000002995",
                    description=None,
                    source=uniprot.GenomeSource.WORMBASE,
                    components=uniprot.ALL_CHROMOSOMES,
                ),
                lineage_info=uniprot.LineageInfo(
                    ncbi_id=6279,
                    species="Brugia malayi",
                    common_name="Filarial nematode worm",
                    tax_string="Eukaryota; Metazoa; Ecdysozoa; Nematoda; Chromadorea; Rhabditida; Spirurina; Spiruromorpha; Filarioidea; Onchocercidae; Brugia.",
                ),
            ),
        ),
        (
            Path("tests/data/UP000006705.xml"),
            uniprot.ProteomeInfo(
                upi="UP000006705",
                taxid=652674,
                is_reference=True,
                is_representative=False,
                proteome_description=None,
                genome_info=uniprot.GenomeInfo(
                    accession="GCA_000861105.1",
                    description="Hepatitis E virus ORF1, ORF2, ORF3,complete cds's; methyltransferase, Y domain, papain-like protease, poly-proline hing",
                    source=uniprot.GenomeSource.ENA,
                    components=uniprot.SelectedComponents.build(
                        [Accession.build("L08816")]
                    ),
                ),
                lineage_info=uniprot.LineageInfo(
                    ncbi_id=652674,
                    species="Hepatitis E virus genotype 1 (isolate Human/China/HeBei/1987)",
                    common_name="HEV",
                    tax_string="Viruses; Riboviria; Orthornavirae; Kitrinoviricota; Alsuviricetes; Hepelivirales; Hepeviridae; Orthohepevirinae; Paslahepevirus; Hepatitis E virus.",
                ),
            ),
        ),
        (
            Path("tests/data/UP000006719.xml"),
            uniprot.ProteomeInfo(
                upi="UP000006719",
                taxid=230407,
                is_reference=True,
                is_representative=True,
                proteome_description=None,
                genome_info=uniprot.GenomeInfo(
                    accession="GCA_000879395.1",
                    description="RefStrain",
                    source=uniprot.GenomeSource.ENA,
                    components=uniprot.SelectedComponents.build(
                        [
                            Accession.build("AB179636"),
                            Accession.build("AB179637"),
                            Accession.build("AB179638"),
                            Accession.build("AB179639"),
                            Accession.build("AB179640"),
                            Accession.build("AB179641"),
                            Accession.build("AB179642"),
                            Accession.build("AB179643"),
                            Accession.build("AY277888"),
                            Accession.build("AY277889"),
                            Accession.build("AY277890"),
                        ]
                    ),
                ),
                lineage_info=uniprot.LineageInfo(
                    ncbi_id=230407,
                    species="Cryphonectria parasitica mycoreovirus 1 (strain 9B21)",
                    common_name="CpMYRV-1",
                    tax_string="Viruses; Riboviria; Orthornavirae; Duplornaviricota; Resentoviricetes; Reovirales; Spinareoviridae; Mycoreovirus; Mycoreovirus 1.",
                ),
            ),
        ),
        (
            Path("tests/data/UP000007078.xml"),
            uniprot.ProteomeInfo(
                upi="UP000007078",
                taxid=766184,
                is_reference=True,
                is_representative=True,
                proteome_description=None,
                genome_info=uniprot.GenomeInfo(
                    accession="AB290918",
                    description="Torque teno midi virus 1 DNA",
                    source=uniprot.GenomeSource.ENA,
                    components=uniprot.ALL_CHROMOSOMES,
                ),
                lineage_info=uniprot.LineageInfo(
                    ncbi_id=766184,
                    species="Torque teno midi virus 1 (isolate MD1-073)",
                    common_name=None,
                    tax_string="Viruses; Anelloviridae; Gammatorquevirus; Torque teno midi virus 1.",
                ),
            ),
        ),
        (
            Path("tests/data/UP000662757.xml"),
            uniprot.ProteomeInfo(
                upi="UP000662757",
                taxid=2813240,
                is_reference=False,
                is_representative=True,
                proteome_description=None,
                genome_info=uniprot.GenomeInfo(
                    accession="GCA_017348405.1",
                    description="Mycobacterium phage prophiGD12-2",
                    source=uniprot.GenomeSource.ENA,
                    components=uniprot.SelectedComponents.build(
                        [
                            Accession.build("MW584207"),
                        ]
                    ),
                ),
                lineage_info=uniprot.LineageInfo(
                    ncbi_id=2813240,
                    species="Mycobacterium phage prophiGD12-2",
                    common_name=None,
                    tax_string="Viruses; Duplodnaviria; Heunggongvirae; Uroviricota; Caudoviricetes; unclassified Caudoviricetes; Siphoviridae; unclassified Siphoviridae.",
                ),
            ),
        ),
        (
            Path("tests/data/UP000000212.xml"),
            uniprot.ProteomeInfo(
                upi="UP000000212",
                taxid=1234679,
                is_reference=False,
                is_representative=True,
                proteome_description=None,
                genome_info=uniprot.GenomeInfo(
                    accession="GCA_000317975.2",
                    description=None,
                    source=uniprot.GenomeSource.ENA,
                    components=uniprot.SelectedComponents.build(
                        [
                            Accession.build("HE999757"),
                        ]
                    ),
                ),
                lineage_info=uniprot.LineageInfo(
                    ncbi_id=1234679,
                    species="Carnobacterium maltaromaticum LMA28",
                    common_name=None,
                    tax_string="Bacteria; Terrabacteria group; Bacillota; Bacilli; Lactobacillales; Carnobacteriaceae; Carnobacterium; Carnobacterium maltaromaticum.",
                ),
            ),
        ),
        (
            Path("tests/data/UP000000224.xml"),
            uniprot.ProteomeInfo(
                upi="UP000000224",
                taxid=1249480,
                is_reference=False,
                is_representative=False,
                proteome_description=None,
                genome_info=uniprot.GenomeInfo(
                    accession="GCA_000310245.1",
                    description=None,
                    source=uniprot.GenomeSource.ENA,
                    components=uniprot.SelectedComponents.build(
                        [
                            Accession.build("CP003921"),
                            Accession.build("CP003920"),
                        ]
                    ),
                ),
                lineage_info=uniprot.LineageInfo(
                    ncbi_id=1249480,
                    species="Candidatus Sulfuricurvum sp. RIFRC-1",
                    common_name=None,
                    tax_string="Bacteria; Pseudomonadota; delta/epsilon subdivisions; Epsilonproteobacteria; Campylobacterales; Thiovulaceae; Sulfuricurvum.",
                ),
            ),
        ),
        (
            Path("tests/data/UP000000230.xml"),
            uniprot.ProteomeInfo(
                upi="UP000000230",
                taxid=399742,
                is_reference=False,
                is_representative=False,
                proteome_description="Enterobacter sp. 638 was isolated in association with poplar, Populus trichocarpa x deltoids, and represents a commonly found endophytic bacterium associated with this tree. Endophytes are bacteria that live within the tissue of a plant without substantively harming it. They can help promote plant growth in several ways, including helping the host overcome toxic effects of environmental pollution. Some Enterobacteriales can affect nodule formation by legumes, fix nitrogen and produce plant hormones. The family Enterobacteriacae includes free-living, commensal and pathogenic bacteria associated with host species ranging from plants to humans (adapted from http://genome.jgi-psf.org/ent_6/ent_6.home.html).",
                genome_info=uniprot.GenomeInfo(
                    accession="GCA_000016325.1",
                    description=None,
                    source=uniprot.GenomeSource.ENA,
                    components=uniprot.SelectedComponents.build(
                        [
                            Accession.build("CP000653"),
                            Accession.build("CP000654"),
                        ]
                    ),
                ),
                lineage_info=uniprot.LineageInfo(
                    ncbi_id=399742,
                    species="Enterobacter sp. (strain 638)",
                    common_name=None,
                    tax_string="Bacteria; Pseudomonadota; Gammaproteobacteria; Enterobacterales; Enterobacteriaceae; Enterobacter; unclassified Enterobacter.",
                ),
            ),
        ),
        (
            Path("tests/data/UP000214863.xml"),
            uniprot.ProteomeInfo(
                upi="UP000214863",
                taxid=2022783,
                is_reference=False,
                is_representative=True,
                proteome_description=None,
                genome_info=uniprot.GenomeInfo(
                    accession="GCA_002210635.1",
                    description="Common bottlenose dolphin gammaherpesvirus 1 strain Sarasota",
                    source=uniprot.GenomeSource.ENA,
                    components=uniprot.SelectedComponents.build(
                        [
                            Accession.build("KY965444"),
                        ]
                    ),
                ),
                lineage_info=uniprot.LineageInfo(
                    ncbi_id=2022783,
                    species="Common bottlenose dolphin gammaherpesvirus 1 strain Sarasota",
                    common_name=None,
                    tax_string="Viruses; Duplodnaviria; Heunggongvirae; Peploviricota; Herviviricetes; Herpesvirales; Herpesviridae; Gammaherpesvirinae; Rhadinovirus; unclassified Rhadinovirus.",
                ),
            ),
        ),
        (
            Path("tests/data/UP000214715.xml"),
            uniprot.ProteomeInfo(
                upi="UP000214715",
                taxid=1497851,
                is_reference=False,
                is_representative=True,
                proteome_description=None,
                genome_info=uniprot.GenomeInfo(
                    accession="GCA_001736955.2",
                    description="Bacillus phage SPG24 DNA, contig00001 sequence.",
                    source=uniprot.GenomeSource.ENA,
                    components=uniprot.SelectedComponents.build(
                        [
                            Accession.build("AB930182"),
                        ]
                    ),
                ),
                lineage_info=uniprot.LineageInfo(
                    ncbi_id=1497851,
                    species="Bacillus phage SPG24",
                    common_name=None,
                    tax_string="Viruses; Duplodnaviria; Heunggongvirae; Uroviricota; Caudoviricetes; Herelleviridae; Bastillevirinae; Nitunavirus; Nitunavirus SPG24.",
                ),
            ),
        ),
        (
            Path("tests/data/UP000214710.xml"),
            uniprot.ProteomeInfo(
                upi="UP000214710",
                taxid=2560515,
                is_reference=False,
                is_representative=True,
                proteome_description=None,
                genome_info=uniprot.GenomeInfo(
                    accession="GCA_002184215.1",
                    description="Grapevine enamovirus-1 isolate SE-BR, partial genome.",
                    source=uniprot.GenomeSource.ENA,
                    components=uniprot.SelectedComponents.build(
                        [
                            Accession.build("KY820716"),
                        ]
                    ),
                ),
                lineage_info=uniprot.LineageInfo(
                    ncbi_id=2560515,
                    species="Grapevine enamovirus 1",
                    common_name=None,
                    tax_string="Viruses; Riboviria; Orthornavirae; Pisuviricota; Pisoniviricetes; Sobelivirales; Solemoviridae; Enamovirus.",
                ),
            ),
        ),
        (
            Path("tests/data/UP000214656.xml"),
            uniprot.ProteomeInfo(
                upi="UP000214656",
                taxid=1654339,
                is_reference=False,
                is_representative=True,
                proteome_description=None,
                genome_info=uniprot.GenomeInfo(
                    accession="GCA_001019975.1",
                    description="Genome",
                    source=uniprot.GenomeSource.ENA,
                    components=uniprot.SelectedComponents.build(
                        [
                            Accession.build("KP774592"),
                            Accession.build("KP774593"),
                        ]
                    ),
                ),
                lineage_info=uniprot.LineageInfo(
                    ncbi_id=1654339,
                    species="Sclerotinia sclerotiorum botybirnavirus 1",
                    common_name=None,
                    tax_string="Viruses; Riboviria; Orthornavirae; Botybirnavirus; unclassified Botybirnavirus.",
                ),
            ),
        ),
        (
            Path("tests/data/UP000214542.xml"),
            uniprot.ProteomeInfo(
                upi="UP000214542",
                taxid=1307954,
                is_reference=False,
                is_representative=False,
                proteome_description=None,
                genome_info=uniprot.GenomeInfo(
                    accession="GCA_000870385.1",
                    description="Maruca vitrata MNPV",
                    source=uniprot.GenomeSource.ENA,
                    components=uniprot.SelectedComponents.build(
                        [
                            Accession.build("EF125867"),
                        ]
                    ),
                ),
                lineage_info=uniprot.LineageInfo(
                    ncbi_id=1307954,
                    species="Maruca vitrata nucleopolyhedrovirus",
                    common_name=None,
                    tax_string="Viruses; Naldaviricetes; Lefavirales; Baculoviridae; Alphabaculovirus.",
                ),
            ),
        ),
        (
            Path("tests/data/UP000001819.xml"),
            uniprot.ProteomeInfo(
                upi="UP000001819",
                taxid=46245,
                is_reference=True,
                is_representative=True,
                proteome_description=None,
                genome_info=uniprot.GenomeInfo(
                    accession="GCF_009870125.1",
                    description=None,
                    source=uniprot.GenomeSource.REF_SEQ,
                    components=uniprot.SelectedComponents.build(
                        [
                            uniprot.UNPLACED,
                            Accession.build("NC_046679"),
                            Accession.build("NC_046680"),
                            Accession.build("NC_046681"),
                            Accession.build("NC_046682"),
                            Accession.build("NC_046683"),
                            Accession.build("NC_046603"),
                        ]
                    ),
                ),
                lineage_info=uniprot.LineageInfo(
                    ncbi_id=46245,
                    species="Drosophila pseudoobscura pseudoobscura",
                    common_name="Fruit fly",
                    tax_string="Eukaryota; Metazoa; Ecdysozoa; Arthropoda; Hexapoda; Insecta; Pterygota; Neoptera; Endopterygota; Diptera; Brachycera; Muscomorpha; Ephydroidea; Drosophilidae; Drosophila; Sophophora.",
                ),
            ),
        ),
        (
            Path("tests/data/UP000007648.xml"),
            uniprot.ProteomeInfo(
                upi="UP000007648",
                taxid=9305,
                is_reference=True,
                is_representative=True,
                proteome_description="The Tasmanian devil (Sarcophilus harrisii) is a member of the family Dasyuridae and is the largest carnivorous marsupial, reaching 76 cm in length and weighing up to 12 kg. It has sharp teeth and strong jaws that can deliver one of the most powerful bites of any mammal. \n\nThe species is restricted to the Australian island of Tasmania and is at risk of extinction in the wild due to an unusual transmissible facial cancer, the devil facial tumor disease (DFTD), which is spread between devils by the transfer of cancer cells on biting.\n\nThe reference proteome of Sarcophilus harrisii is derived from the genome sequence prepared by the Tasmanian Devil Genome Project by the Center for Comparative Genomics and Bioinformatics at Pennsylvania State University. The project analyzed genomes of two Tasmanian devil individuals, one healthy animal and one with DFTD. The size of the genome is about 3.3 Gb.",
                genome_info=uniprot.GenomeInfo(
                    accession="GCA_902635505.1",
                    description=None,
                    components=uniprot.SelectedComponents.build(
                        [
                            WgsPrefix.build("AFEY01000000"),
                        ]
                    ),
                    source=uniprot.GenomeSource.ENSEMBL,
                ),
                lineage_info=uniprot.LineageInfo(
                    ncbi_id=9305,
                    species="Sarcophilus harrisii",
                    tax_string="Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Metatheria; Dasyuromorphia; Dasyuridae; Sarcophilus.",
                    common_name="Tasmanian devil",
                ),
            ),
        ),
        (
            Path("tests/data/UP000019116.xml"),
            uniprot.ProteomeInfo(
                upi="UP000019116",
                taxid=4565,
                is_reference=False,
                is_representative=True,
                proteome_description=None,
                genome_info=uniprot.GenomeInfo(
                    accession="GCA_900519105.1",
                    description=None,
                    source=uniprot.GenomeSource.ENSEMBL_PLANTS,
                    components=uniprot.SelectedComponents.build(
                        [
                            Accession.build("AB042240"),
                            Accession.build("AP008982"),
                            uniprot.UNPLACED,
                        ]
                    ),
                ),
                lineage_info=uniprot.LineageInfo(
                    ncbi_id=4565,
                    species="Triticum aestivum",
                    common_name="Wheat",
                    tax_string="Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta; Spermatophyta; Magnoliopsida; Liliopsida; Poales; Poaceae; BOP clade; Pooideae; Triticodae; Triticeae; Triticinae; Triticum.",
                ),
            ),
        ),
        (
            Path("tests/data/UP000093657.xml"),
            uniprot.ProteomeInfo(
                upi="UP000093657",
                taxid=1849840,
                is_reference=False,
                is_representative=True,
                proteome_description=None,
                genome_info=uniprot.GenomeInfo(
                    accession="GCA_001675455.1",
                    description=None,
                    source=uniprot.GenomeSource.ENA,
                    components=uniprot.SelectedComponents.build(
                        [
                            WgsPrefix.build("LXWX01000000"),
                        ]
                    ),
                ),
                lineage_info=uniprot.LineageInfo(
                    ncbi_id=1849840,
                    species="Methylosinus sp. 3S-1",
                    common_name=None,
                    tax_string="Bacteria; Pseudomonadota; Alphaproteobacteria; Hyphomicrobiales; Methylocystaceae; Methylosinus; unclassified Methylosinus.",
                ),
            ),
        ),
        (
            Path("tests/data/UP000214540.xml"),
            uniprot.ProteomeInfo(
                upi="UP000214540",
                taxid=1983562,
                is_reference=False,
                is_representative=True,
                proteome_description=None,
                genome_info=uniprot.GenomeInfo(
                    accession="GCA_002237195.1",
                    description="Lake Sinai Virus SA1 ORF1, RNA-dependent RNA polymerase, and ORF4 genes",
                    source=uniprot.GenomeSource.ENA,
                    components=uniprot.SelectedComponents.build(
                        [
                            Accession.build("KY354243"),
                        ]
                    ),
                ),
                lineage_info=uniprot.LineageInfo(
                    ncbi_id=1983562,
                    species="Lake Sinai Virus SA1",
                    common_name=None,
                    tax_string="Viruses; Riboviria; Orthornavirae; Kitrinoviricota; Magsaviricetes; Nodamuvirales; Sinhaliviridae; Sinaivirus; unclassified Sinaivirus.",
                ),
            ),
        ),
        (
            Path("tests/data/UP000035680.xml"),
            uniprot.ProteomeInfo(
                upi="UP000035680",
                taxid=75913,
                is_reference=False,
                is_representative=True,
                proteome_description=None,
                genome_info=uniprot.GenomeInfo(
                    accession="GCA_001028725.1",
                    description=None,
                    source=uniprot.GenomeSource.WORMBASE,
                    components=uniprot.ALL_CHROMOSOMES,
                ),
                lineage_info=uniprot.LineageInfo(
                    ncbi_id=75913,
                    species="Strongyloides venezuelensis",
                    common_name="Threadworm",
                    tax_string="Eukaryota; Metazoa; Ecdysozoa; Nematoda; Chromadorea; Rhabditida; Tylenchina; Panagrolaimomorpha; Strongyloidoidea; Strongyloididae; Strongyloides.",
                ),
            ),
        ),
        (
            Path("tests/data/UP000000408.xml"),
            uniprot.ProteomeInfo(
                upi="UP000000408",
                taxid=650136,
                is_reference=True,
                is_representative=False,
                proteome_description=None,
                genome_info=uniprot.GenomeInfo(
                    accession="GCA_000857805.1",
                    description="Porcine enterovirus 1 serotype 1 RNA for complete polyprotein, isolate F65",
                    source=uniprot.GenomeSource.ENA,
                    components=uniprot.SelectedComponents.build(
                        [Accession.build("AJ011380")]
                    ),
                ),
                lineage_info=uniprot.LineageInfo(
                    ncbi_id=650136,
                    species="Porcine teschovirus 1 (isolate Pig/United Kingdom/F65/1967)",
                    common_name="PTV-1",
                    tax_string="Viruses; Riboviria; Orthornavirae; Pisuviricota; Pisoniviricetes; Picornavirales; Picornaviridae; Caphthovirinae; Teschovirus; Teschovirus A; teschovirus A1.",
                ),
            ),
        ),
        (
            Path("tests/data/UP000008177.xml"),
            uniprot.ProteomeInfo(
                upi="UP000008177",
                taxid=999810,
                is_reference=False,
                is_representative=True,
                proteome_description=None,
                genome_info=uniprot.GenomeInfo(
                    accession="GCA_000292645.1",
                    description=None,
                    source=uniprot.GenomeSource.ENA,
                    components=uniprot.SelectedComponents.build(
                        [
                            uniprot.UNPLACED,
                            Accession.build("FQ790245"),
                            Accession.build("FQ790246"),
                            Accession.build("FQ790247"),
                            Accession.build("FQ790248"),
                            Accession.build("FQ790249"),
                            Accession.build("FQ790250"),
                            Accession.build("FQ790251"),
                            Accession.build("FQ790252"),
                            Accession.build("FQ790253"),
                            Accession.build("FQ790254"),
                            Accession.build("FQ790255"),
                            Accession.build("FQ790256"),
                            Accession.build("FQ790257"),
                            Accession.build("FQ790258"),
                            Accession.build("FQ790259"),
                            Accession.build("FQ790260"),
                            Accession.build("FQ790261"),
                            Accession.build("FQ790262"),
                            Accession.build("FQ790263"),
                            Accession.build("FQ790264"),
                            Accession.build("FQ790265"),
                            Accession.build("FQ790266"),
                            Accession.build("FQ790267"),
                            Accession.build("FQ790268"),
                            Accession.build("FQ790269"),
                            Accession.build("FQ790270"),
                            Accession.build("FQ790271"),
                            Accession.build("FQ790272"),
                            Accession.build("FQ790273"),
                            Accession.build("FQ790274"),
                            Accession.build("FQ790275"),
                            Accession.build("FQ790276"),
                            Accession.build("FQ790277"),
                            Accession.build("FQ790278"),
                            Accession.build("FQ790279"),
                            Accession.build("FQ790280"),
                            Accession.build("FQ790281"),
                            Accession.build("FQ790282"),
                            Accession.build("FQ790283"),
                            Accession.build("FQ790284"),
                            Accession.build("FQ790285"),
                            Accession.build("FQ790286"),
                            Accession.build("FQ790287"),
                            Accession.build("FQ790288"),
                            Accession.build("FQ790289"),
                            Accession.build("FQ790290"),
                            Accession.build("FQ790291"),
                            Accession.build("FQ790292"),
                            Accession.build("FQ790293"),
                            Accession.build("FQ790294"),
                            Accession.build("FQ790295"),
                            Accession.build("FQ790296"),
                            Accession.build("FQ790297"),
                            Accession.build("FQ790298"),
                            Accession.build("FQ790299"),
                            Accession.build("FQ790300"),
                            Accession.build("FQ790301"),
                            Accession.build("FQ790302"),
                            Accession.build("FQ790303"),
                            Accession.build("FQ790304"),
                            Accession.build("FQ790305"),
                            Accession.build("FQ790306"),
                            Accession.build("FQ790307"),
                            Accession.build("FQ790308"),
                            Accession.build("FQ790309"),
                            Accession.build("FQ790310"),
                            Accession.build("FQ790311"),
                            Accession.build("FQ790312"),
                            Accession.build("FQ790313"),
                            Accession.build("FQ790314"),
                            Accession.build("FQ790315"),
                            Accession.build("FQ790316"),
                            Accession.build("FQ790317"),
                            Accession.build("FQ790318"),
                            Accession.build("FQ790319"),
                            Accession.build("FQ790320"),
                            Accession.build("FQ790321"),
                            Accession.build("FQ790322"),
                            Accession.build("FQ790323"),
                            Accession.build("FQ790324"),
                            Accession.build("FQ790325"),
                            Accession.build("FQ790326"),
                            Accession.build("FQ790327"),
                            Accession.build("FQ790328"),
                            Accession.build("FQ790329"),
                            Accession.build("FQ790330"),
                            Accession.build("FQ790331"),
                            Accession.build("FQ790332"),
                            Accession.build("FQ790333"),
                            Accession.build("FQ790334"),
                            Accession.build("FQ790335"),
                            Accession.build("FQ790336"),
                            Accession.build("FQ790337"),
                            Accession.build("FQ790338"),
                            Accession.build("FQ790339"),
                            Accession.build("FQ790340"),
                            Accession.build("FQ790341"),
                            Accession.build("FQ790342"),
                            Accession.build("FQ790343"),
                            Accession.build("FQ790344"),
                            Accession.build("FQ790345"),
                            Accession.build("FQ790346"),
                            Accession.build("FQ790347"),
                            Accession.build("FQ790348"),
                            Accession.build("FQ790349"),
                            Accession.build("FQ790350"),
                            Accession.build("FQ790351"),
                            Accession.build("FQ790352"),
                            Accession.build("FQ790353"),
                            Accession.build("FQ790354"),
                            Accession.build("FQ790355"),
                            Accession.build("FQ790356"),
                            Accession.build("FQ790357"),
                            Accession.build("FQ790358"),
                            Accession.build("FQ790359"),
                            Accession.build("FQ790360"),
                            Accession.build("FQ790361"),
                            Accession.build("FQ790362"),
                            WgsPrefix.build("ALOC01"),
                        ]
                    ),
                ),
                lineage_info=uniprot.LineageInfo(
                    ncbi_id=999810,
                    species="Botryotinia fuckeliana (strain T4)",
                    common_name="Noble rot fungus",
                    tax_string="Eukaryota; Fungi; Dikarya; Ascomycota; Pezizomycotina; Leotiomycetes; Helotiales; Sclerotiniaceae; Botrytis.",
                ),
            ),
        ),
        (
            Path("tests/data/UP000035681.xml"),
            uniprot.ProteomeInfo(
                upi="UP000035681",
                taxid=6248,
                is_reference=False,
                is_representative=True,
                proteome_description=None,
                genome_info=uniprot.GenomeInfo(
                    accession="GCA_000947215.1",
                    description=None,
                    components=uniprot.ALL_CHROMOSOMES,
                    source=uniprot.GenomeSource.WORMBASE,
                ),
                lineage_info=uniprot.LineageInfo(
                    ncbi_id=6248,
                    species="Strongyloides stercoralis",
                    common_name="Threadworm",
                    tax_string="Eukaryota; Metazoa; Ecdysozoa; Nematoda; Chromadorea; Rhabditida; Tylenchina; Panagrolaimomorpha; Strongyloidoidea; Strongyloididae; Strongyloides.",
                ),
            ),
        ),
        (
            Path("test/data/UP000006540.json"),
            None,  # TODO This is a nice example of a null accesion to handle
        ),
    ],
)
def test_parses_expected_data(path: Path, expected: uniprot.ProteomeInfo):
    data = list(uniprot.proteomes(path, set()))
    assert len(data) == 1
    assert data[0] == expected


@pytest.mark.skip()
@pytest.mark.parametrize(
    "path,ignore,expected",
    [
        (Path("tests/data/proteome.xml"), set(), 21739),
        (Path("tests/data/proteome.xml"), {"UP000005640"}, 21738),
    ],
)
def test_parses_expected_number(path: Path, ignore: ty.Set[str], expected: int):
    assert len(list(uniprot.proteomes(path, ignore))) == expected


@pytest.mark.parametrize(
    "accessions,given,expected",
    [
        (
            [Accession.build("J01415"), uniprot.UNPLACED],
            Accession.build("CM000663"),
            [],
        ),
        (
            [Accession.build("J01415"), uniprot.UNPLACED],
            Accession.build("J01415"),
            [Accession.build("J01415")],
        ),
        (
            [Accession.build("J01415", aliases=(Accession.build("CM000663"),))],
            Accession.build("CM000663"),
            [Accession.build("J01415", aliases=(Accession.build("CM000663"),))],
        ),
        (
            [Accession.build("J01415", aliases=(Accession.build("CM000663"),))],
            Accession.build("J01416"),
            [],
        ),
        (
            [Accession.build("J01415", aliases=(Accession.build("CM000663"),))],
            Accession.build("CM000664"),
            [],
        ),
        (
            [Accession.build("CM000663")],
            Accession.build("J01415", aliases=(Accession.build("CM000663"),)),
            [Accession.build("CM000663")],
        ),
        ([uniprot.UNPLACED], Accession.build("J01415"), []),
        ([uniprot.UNPLACED], uniprot.UNPLACED, [uniprot.UNPLACED]),
        ([Accession.build("J01415")], uniprot.UNPLACED, []),
        (
            [WgsPrefix.build("ANNZ01000000"), WgsSequenceId.build("ANNZ01000001")],
            Accession.build("CM000663"),
            [],
        ),
        (
            [WgsPrefix.build("ANNZ01000000"), WgsSequenceId.build("ANNZ01000001")],
            WgsPrefix.build("ANNZ01"),
            [WgsPrefix.build("ANNZ01000000"), WgsSequenceId.build("ANNZ01000001")],
        ),
        (
            [WgsPrefix.build("ANNZ01000000"), WgsSequenceId.build("ANNZ01000001")],
            WgsSequenceId.build("ANNZ01000001"),
            [WgsPrefix.build("ANNZ01000000"), WgsSequenceId.build("ANNZ01000001")],
        ),
        (
            [WgsPrefix.build("ANNZ01000000"), WgsSequenceId.build("ANNZ01000001")],
            WgsSequenceId.build("ANNZ01000002"),
            [WgsPrefix.build("ANNZ01000000")],
        ),
        (
            [WgsPrefix.build("ANNZ01000000"), WgsSequenceId.build("ANNZ01000001")],
            WgsSequenceId.build("ANNZ02000001"),
            [WgsPrefix.build("ANNZ01000000"), WgsSequenceId.build("ANNZ01000001")],
        ),
    ],
)
def test_selected_components_can_get_expected_matching_components(
    accessions, given, expected
):
    selected = uniprot.SelectedComponents.build(accessions)
    assert selected.matching(given) == expected
