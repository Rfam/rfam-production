"""
Copyright [2009-2016] EMBL-European Bioinformatics Institute
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

import os
import shutil

from scripts.export.genomes import genome_fetch as gf


# --------------------------------------------------------------------------------------------------
def test_find_proteomes_without_accessions():
    id_pairs = gf.load_upid_gca_pairs()

    for upid in id_pairs.keys():
        print "%s %s" % (upid, id_pairs[upid])

        """
        accs = fetch_genome_accessions(upid, id_pairs[upid])

        if len(accs) == 0:
            print upid
        else:
            print "%s %s\n"%(upid,str(accs))
        """


# --------------------------------------------------------------------------------------------------

def test_fetch_ref_proteomes():
    ref_prots = None
    ref_prots = gf.fetch_ref_proteomes()

    assert len(ref_prots) != 0 or ref_prots is not None


# --------------------------------------------------------------------------------------------------

def test_export_gca_accessions():
    upid_gca_pairs = None
    upid_gca_pairs = gf.export_gca_accessions("input/UPID_GCA.tsv")

    assert upid_gca_pairs is not None
    assert len(upid_gca_pairs.keys()) != 0


# --------------------------------------------------------------------------------------------------
def test_extract_genome_acc():
    # homo sapiens
    gen_acc = None
    gen_acc = gf.extract_genome_acc("http://www.uniprot.org/proteomes/UP000005640.rdf")

    print gen_acc
    assert gen_acc == -1 or gen_acc is not None


# --------------------------------------------------------------------------------------------------
def test_proteome_rdf_scanner():
    # homo sapiens
    accs = None
    accs = gf.proteome_rdf_scanner("http://www.uniprot.org/proteomes/UP000005640.rdf")

    print accs
    assert len(accs) != 0 or accs is not None


# --------------------------------------------------------------------------------------------------
def test_fetch_genome_acc():
    gen_accs = None
    gen_accs = gf.fetch_genome_acc("UP000005640")

    print gen_accs
    assert gen_accs is not None or len(gen_accs) != 0


# --------------------------------------------------------------------------------------------------
def test_fetch_ena_file():
    # need something smaller here
    if not os.path.exists("/tmp/gen_test"):
        os.mkdir("/tmp/gen_test")

    # Citrus psorosis virus
    check = gf.fetch_ena_file("AY654894", "fasta", "/tmp/gen_test")

    assert check is True
    shutil.rmtree("/tmp/gen_test")


# --------------------------------------------------------------------------------------------------
def test_extract_assembly_accs():
    gen_accs = None
    gen_accs = gf.extract_assembly_accs("GCA_000001405.23")

    print gen_accs
    assert gen_accs is not None or len(gen_accs) != 0


# --------------------------------------------------------------------------------------------------
"""
def test_fetch_genome():
    gen = "GCA_000320365.1"
    dest_dir = "/tmp/gen_test"

    gf.fetch_genome(gen, dest_dir)
    gen_files = os.listdir(os.path.join("/tmp/gen_test", gen.partition('.')[0]))

    assert len(gen_files) != 0
    shutil.rmtree(os.path.join("/tmp/gen_test", gen.partition('.')[0]))
"""


# --------------------------------------------------------------------------------------------------
def test_rdf_accession_search():
    # homo sapiens
    rdf_accs = None
    rdf_accs = gf.rdf_accession_search("UP000005640", "/embl/")

    assert rdf_accs is not None or len(rdf_accs) != 0


# --------------------------------------------------------------------------------------------------
def test_load_upid_gca_file():
    id_pairs = None
    id_pairs = gf.load_upid_gca_file("input/UPID_GCA.tsv")

    print id_pairs
    assert id_pairs is not None or len(id_pairs.keys()) != 0


# --------------------------------------------------------------------------------------------------
def test_load_upid_gca_pairs():
    id_pairs = None
    id_pairs = gf.load_upid_gca_pairs()

    print id_pairs
    assert id_pairs is not None or len(id_pairs.keys()) != 0


# --------------------------------------------------------------------------------------------------
def test_fetch_genome_accessions():
    gen_accs = None
    # need to get all the cases here - Unittests test_cases needed
    gen_accs = gf.fetch_genome_accessions("UP000005640", "GCA_000001405.23")

    print gen_accs
    assert gen_accs is not None or len(gen_accs) != 0


# --------------------------------------------------------------------------------------------------

def test_assembly_report_parser():
    accessions = None
    # homo sapiens
    report_url = "ftp://ftp.ebi.ac.uk/pub/databases/ena/assembly/GCA_000/GCA_000001/GCA_000001405.23_sequence_report.txt"
    accessions = gf.assembly_report_parser(report_url)

    print accessions
    assert len(accessions) != 0 or accessions is not None


# --------------------------------------------------------------------------------------------------
def test_get_wgs_set_accession():
    wgs_acc = None
    wgs_acc = gf.get_wgs_set_accession("AYNF", "1")

    print wgs_acc
    assert wgs_acc is not None


# --------------------------------------------------------------------------------------------------
def test_fetch_wgs_range_accs():
    wgs_accs = None
    # Escherichia coli LAU-EC6
    wgs_accs = gf.fetch_wgs_range_accs("AYNF01000001-AYNF01000106")

    print wgs_accs
    assert wgs_accs is not None or len(wgs_accs) != 0


# --------------------------------------------------------------------------------------------------

def test_download_genomes():
    # Moraxella macacae
    gen = "GCA_000320365.1"
    dest_dir = "/tmp/gen_test"

    if not os.path.exists(dest_dir):
        os.mkdir(dest_dir)

    check = gf.download_genomes(gen, dest_dir)

    gen_files = os.listdir(os.path.join(dest_dir, gen.partition('.')[0]))

    print gen_files
    assert check is not None and len(gen_files) != 0

    shutil.rmtree(dest_dir)


# --------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    # test_find_proteomes_without_accessions()
    # test_fetch_ref_proteomes()
    # test_extract_genome_acc()
    # test_proteome_rdf_scanner()
    # test_fetch_genome_acc()
    # test_extract_assembly_accs()
    # test_fetch_wgs_range_accs()
    # test_fetch_genome_accessions()
    # test_load_upid_gca_pairs()
    # test_load_upid_gca_file()
    # test_get_wgs_set_accession()
    # test_assembly_report_parser()
    # test_download_genomes()
    test_fetch_genome()
    pass
