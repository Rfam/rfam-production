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
    "path,expected",
    [
        (
            Path("tests/data/UP000005640.xml"),
            uniprot.ProteomeInfo(
                upi="UP000005640",
                taxid="9606",
                is_reference=True,
                is_representative=True,
                proteome_description='''Homo sapiens (Homo sapiens sapiens) or modern humans are the only living species of the evolutionary branch of great apes known as hominids. Divergence of early humans from chimpanzees and gorillas is estimated to have occurred between 4 and 8 million years ago. The genus Homo (Homo habilis) appeared in Africa around 2.3 million years ago and shows the first signs of stone tool usage. The exact lineage of Homo species ie: H. habilis/H. ergaster to H. erectus to H. rhodesiensis/H.heidelbergensis to H. sapiens is still hotly disputed. However, continuing evolution and in particular larger brain size and complexity culminates in Homo sapiens. The first anatomically modern humans appear in the fossil record around 200,000 years ago. Modern humans migrated across the globe essentially as hunter-gatherers until around 12,000 years ago when the practice of agriculture and animal domestication enabled large populations to grow leading to the development of civilizations.''',
                genome_info=uniprot.GenomeInfo(
                    accession="GCA_000001405.27",
                    description=None,
                    components=uniprot.SelectedComponents(
                        accessions=[
                            "CM000663",
                            "CM000664",
                            "CM000665",
                            "CM000666",
                            "CM000667",
                            "CM000668",
                            "CM000669",
                            "CM000670",
                            "CM000671",
                            "CM000672",
                            "CM000673",
                            "CM000674",
                            "CM000675",
                            "CM000676",
                            "CM000677",
                            "CM000678",
                            "CM000679",
                            "CM000680",
                            "CM000681",
                            "CM000682",
                            "CM000683",
                            "CM000684",
                            "CM000685",
                            "CM000686",
                            "J01415",
                            uniprot.UNPLACED,
                        ]
                    ),
                ),
            ),
        ),
        (
            Path("tests/data/UP000000226.xml"),
            uniprot.ProteomeInfo(
                upi="UP000000226",
                taxid="3885",
                is_reference=False,
                is_representative=True,
                proteome_description='The common bean is most widely cultivated of all beans in temperate regions, and is widely cultivated in semitropical regions. In some regions of the world dry beans furnish a large portion of the protein needs of low and middle class families.',
                genome_info=uniprot.GenomeInfo(
                    accession="GCA_000499845.1",
                    description=None,
                    components=uniprot.SelectedComponents(
                        accessions=[
                            "CM002288",
                            "CM002289",
                            "CM002290",
                            "CM002291",
                            "CM002292",
                            "CM002293",
                            "CM002294",
                            "CM002295",
                            "CM002296",
                            "CM002297",
                            "CM002298",
                            "ANNZ01000000",
                        ]
                    ),
                ),
            ),
        ),
        (
            Path("tests/data/UP000006672.xml"),
            uniprot.ProteomeInfo(
                upi="UP000006672",
                taxid="6279",
                is_reference=True,
                is_representative=True,
                proteome_description=None,
                genome_info=uniprot.GenomeInfo(
                    accession="GCA_000002995",
                    description=None,
                    components=uniprot.ALL_CHROMOSOMES,
                ),
            ),
        ),
        (
            Path("tests/data/UP000006705.xml"),
            uniprot.ProteomeInfo(
                upi="UP000006705",
                taxid="652674",
                is_reference=True,
                is_representative=False,
                proteome_description=None,
                genome_info=uniprot.GenomeInfo(
                    accession="GCA_000861105.1",
                    description="Hepatitis E virus ORF1, ORF2, ORF3,complete cds's; methyltransferase, Y domain, papain-like protease, poly-proline hing",
                    components=uniprot.SelectedComponents(accessions=["L08816"]),
                ),
            ),
        ),
        (
            Path("tests/data/UP000006719.xml"),
            uniprot.ProteomeInfo(
                upi="UP000006719",
                taxid="230407",
                is_reference=True,
                is_representative=True,
                proteome_description=None,
                genome_info=uniprot.GenomeInfo(
                    accession="GCA_000879395.1",
                    description='RefStrain',
                    components=uniprot.SelectedComponents(
                        accessions=[
                            "AB179636",
                            "AB179637",
                            "AB179638",
                            "AB179639",
                            "AB179640",
                            "AB179641",
                            "AB179642",
                            "AB179643",
                            "AY277888",
                            "AY277889",
                            "AY277890",
                        ]
                    ),
                ),
            ),
        ),
        (
            Path("tests/data/UP000007078.xml"),
            uniprot.ProteomeInfo(
                upi="UP000007078",
                taxid="766184",
                is_reference=True,
                is_representative=True,
                proteome_description=None,
                genome_info=uniprot.GenomeInfo(
                    accession="AB290918",
                    description='Torque teno midi virus 1 DNA',
                    components=uniprot.ALL_CHROMOSOMES,
                ),
            ),
        ),
        (
            Path("tests/data/UP000662757.xml"),
            uniprot.ProteomeInfo(
                upi="UP000662757",
                taxid="2813240",
                is_reference=False,
                is_representative=True,
                proteome_description=None,
                genome_info=uniprot.GenomeInfo(
                    accession="GCA_017348405.1",
                    description='Mycobacterium phage prophiGD12-2',
                    components=uniprot.SelectedComponents(
                        accessions=[
                            "MW584207",
                        ]
                    ),
                ),
            ),
        ),
        (
            Path("tests/data/UP000000212.xml"),
            uniprot.ProteomeInfo(
                upi="UP000000212",
                taxid="1234679",
                is_reference=False,
                is_representative=True,
                proteome_description=None,
                genome_info=uniprot.GenomeInfo(
                    accession="GCA_000317975.2",
                    description='Carnobacterium maltaromaticum LMA28 complete genome',
                    components=uniprot.SelectedComponents(
                        accessions=[
                            "HE999757",
                        ]
                    ),
                ),
            ),
        ),
        (
            Path("tests/data/UP000000224.xml"),
            uniprot.ProteomeInfo(
                upi="UP000000224",
                taxid="1249480",
                is_reference=False,
                is_representative=False,
                proteome_description=None,
                genome_info=uniprot.GenomeInfo(
                    accession="GCA_000310245.1",
                    description=None,
                    components=uniprot.SelectedComponents(
                        accessions=[
                            "CP003921",
                            "CP003920",
                        ]
                    ),
                ),
            ),
        ),
        (
            Path("tests/data/UP000000230.xml"),
            uniprot.ProteomeInfo(
                upi="UP000000230",
                taxid="399742",
                is_reference=False,
                is_representative=False,
                proteome_description='Enterobacter sp. 638 was isolated in association with poplar, Populus trichocarpa x deltoids, and represents a commonly found endophytic bacterium associated with this tree. Endophytes are bacteria that live within the tissue of a plant without substantively harming it. They can help promote plant growth in several ways, including helping the host overcome toxic effects of environmental pollution. Some Enterobacteriales can affect nodule formation by legumes, fix nitrogen and produce plant hormones. The family Enterobacteriacae includes free-living, commensal and pathogenic bacteria associated with host species ranging from plants to humans (adapted from http://genome.jgi-psf.org/ent_6/ent_6.home.html).',
                genome_info=uniprot.GenomeInfo(
                    accession="GCA_000016325.1",
                    description=None,
                    components=uniprot.SelectedComponents(
                        accessions=[
                            "CP000653",
                            "CP000654",
                        ]
                    ),
                ),
            ),
        ),
        (
            Path("tests/data/UP000214863.xml"),
            uniprot.ProteomeInfo(
                upi="UP000214863",
                taxid="2022783",
                is_reference=False,
                is_representative=True,
                proteome_description=None,
                genome_info=uniprot.GenomeInfo(
                    accession="GCA_002210635.1",
                    description='Common bottlenose dolphin gammaherpesvirus 1 strain Sarasota',
                    components=uniprot.SelectedComponents(
                        accessions=[
                            "KY965444",
                        ]
                    ),
                ),
            ),
        ),
        (
            Path("tests/data/UP000214715.xml"),
            uniprot.ProteomeInfo(
                upi="UP000214715",
                taxid="1497851",
                is_reference=False,
                is_representative=True,
                proteome_description=None,
                genome_info=uniprot.GenomeInfo(
                    accession="GCA_001736955.2",
                    description='Bacillus phage SPG24 DNA, contig00001 sequence.',
                    components=uniprot.SelectedComponents(
                        accessions=[
                            "AB930182",
                        ]
                    ),
                ),
            ),
        ),
        (
            Path("tests/data/UP000214710.xml"),
            uniprot.ProteomeInfo(
                upi="UP000214710",
                taxid="2560515",
                is_reference=False,
                is_representative=True,
                proteome_description=None,
                genome_info=uniprot.GenomeInfo(
                    accession="GCA_002184215.1",
                    description='Grapevine enamovirus-1 isolate SE-BR, partial genome.',
                    components=uniprot.SelectedComponents(
                        accessions=[
                            "KY820716",
                        ]
                    ),
                ),
            ),
        ),
        (
            Path("tests/data/UP000214656.xml"),
            uniprot.ProteomeInfo(
                upi="UP000214656",
                taxid="1654339",
                is_reference=False,
                is_representative=True,
                proteome_description=None,
                genome_info=uniprot.GenomeInfo(
                    accession="GCA_001019975.1",
                    description='Genome',
                    components=uniprot.SelectedComponents(
                        accessions=[
                            "KP774592",
                            "KP774593",
                        ]
                    ),
                ),
            ),
        ),
        (
            Path("tests/data/UP000214542.xml"),
            uniprot.ProteomeInfo(
                upi="UP000214542",
                taxid="1307954",
                is_reference=False,
                is_representative=False,
                proteome_description=None,
                genome_info=uniprot.GenomeInfo(
                    accession="GCA_000870385.1",
                    description='Maruca vitrata MNPV',
                    components=uniprot.SelectedComponents(
                        accessions=[
                            "EF125867",
                        ]
                    ),
                ),
            ),
        ),
        (
            Path("tests/data/UP000001819.xml"),
            uniprot.ProteomeInfo(
                upi="UP000001819",
                taxid="46245",
                is_reference=True,
                is_representative=True,
                proteome_description=None,
                genome_info=uniprot.GenomeInfo(
                    accession="GCF_009870125.1",
                    description=None,
                    components=uniprot.SelectedComponents(
                        accessions=[
                            uniprot.UNPLACED,
                            "NC_046679",
                            "NC_046680",
                            "NC_046681",
                            "NC_046682",
                            "NC_046683",
                            "NC_046603",
                        ]
                    ),
                ),
            ),
        ),
        (
            Path("tests/data/UP000007648.xml"),
            uniprot.ProteomeInfo(
                upi="UP000007648",
                taxid="9305",
                is_reference=True,
                is_representative=True,
                proteome_description='''The Tasmanian devil (Sarcophilus harrisii) is a member of the family Dasyuridae and is the largest carnivorous marsupial, reaching 76 cm in length and weighing up to 12 kg. It has sharp teeth and strong jaws that can deliver one of the most powerful bites of any mammal. 

The species is restricted to the Australian island of Tasmania and is at risk of extinction in the wild due to an unusual transmissible facial cancer, the devil facial tumor disease (DFTD), which is spread between devils by the transfer of cancer cells on biting.

The reference proteome of Sarcophilus harrisii is derived from the genome sequence prepared by the Tasmanian Devil Genome Project by the Center for Comparative Genomics and Bioinformatics at Pennsylvania State University. The project analyzed genomes of two Tasmanian devil individuals, one healthy animal and one with DFTD. The size of the genome is about 3.3 Gb.''',
                genome_info=uniprot.GenomeInfo(
                    accession="GCA_902635505.1",
                    description=None,
                    components=uniprot.SelectedComponents(
                        accessions=[
                            "AFEY01000000",
                        ]
                    ),
                ),
            ),
        ),
        (
            Path("tests/data/UP000019116.xml"),
            uniprot.ProteomeInfo(
                upi="UP000019116",
                taxid="4565",
                is_reference=False,
                is_representative=True,
                proteome_description=None,
                genome_info=uniprot.GenomeInfo(
                    accession="GCA_900519105.1",
                    description=None,
                    components=uniprot.SelectedComponents(
                        accessions=[
                            "AB042240",
                            "AP008982",
                            uniprot.UNPLACED,
                        ]
                    ),
                ),
            ),
        ),
        (
            Path("tests/data/UP000093657.xml"),
            uniprot.ProteomeInfo(
                upi="UP000093657",
                taxid="1849840",
                is_reference=False,
                is_representative=True,
                proteome_description=None,
                genome_info=uniprot.GenomeInfo(
                    accession="GCA_001675455.1",
                    description=None,
                    components=uniprot.SelectedComponents(
                        accessions=[
                            "LXWX01000000",
                        ]
                    ),
                ),
            ),
        ),
        (
            Path("tests/data/UP000214540.xml"),
            uniprot.ProteomeInfo(
                upi="UP000214540",
                taxid="1983562",
                is_reference=False,
                is_representative=True,
                proteome_description=None,
                genome_info=uniprot.GenomeInfo(
                    accession="GCA_002237195.1",
                    description='Lake Sinai Virus SA1 ORF1, RNA-dependent RNA polymerase, and ORF4 genes',
                    components=uniprot.SelectedComponents(
                        accessions=[
                            "KY354243",
                        ]
                    ),
                ),
            ),
        ),
        (
            Path("tests/data/UP000035680.xml"),
            uniprot.ProteomeInfo(
                upi="UP000035680",
                taxid="75913",
                is_reference=False,
                is_representative=True,
                proteome_description=None,
                genome_info=uniprot.GenomeInfo(
                    accession="GCA_001028725.1",
                    description=None,
                    components=uniprot.ALL_CHROMOSOMES,
                ),
            ),
        ),
        (
            Path("tests/data/UP000000408.xml"),
            uniprot.ProteomeInfo(
                upi="UP000000408",
                taxid="650136",
                is_reference=True,
                is_representative=False,
                proteome_description=None,
                genome_info=uniprot.GenomeInfo(
                    accession="GCA_000857805.1",
                    description='Porcine enterovirus 1 serotype 1 RNA for complete polyprotein, isolate F65',
                    components=uniprot.SelectedComponents(accessions=['AJ011380'])
                ),
            ),
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
