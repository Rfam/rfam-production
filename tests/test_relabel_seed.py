import unittest
import scripts.preprocessing.relabel_seed as rs

# -----------------------------------------------------------------------


class TestRelabelSeedFunctions(unittest.TestCase):

    def test_load_fasta_file_to_dict(self):
        pass

    # -----------------------------------------------------------------------

    def test_sequence_to_md5(self):

        "https://rnacentral.org/api/v1/rna/URS0000677DD8.json"

        # RNAcentral accession: URS0000677DD8
        sequence = "UCAAUAAUGAAAUCUUCUGAUUUGGUGAGAAAUAAUGCCUUAAAAUUACACUCAAUAGGAUUAUGCUGAGG"

        md5_hash = rs.sequence_to_md5(sequence)

        self.assertEqual(md5_hash, "e25b1f5967d21e39e394edb4d3054d0e")

    # -----------------------------------------------------------------------

    def test_fetch_RNAcentral_id(self):

        sequence = "UCAAUAAUGAAAUCUUCUGAUUUGGUGAGAAAUAAUGCCUUAAAAUUACACUCAAUAGGAUUAUGCUGAGG"

        rnacentral_id = rs.fetch_RNAcentral_id(sequence)

        self.assertEqual(rnacentral_id, "URS0000677DD8")

    # -----------------------------------------------------------------------

    def test_generate_RNAcentral_seed_id(self):

        sequence = "UCAAUAAUGAAAUCUUCUGAUUUGGUGAGAAAUAAUGCCUUAAAAUUACACUCAAUAGGAUUAUGCUGAGG"

        seed_id = rs.generate_RNAcentral_seed_id(sequence)

        self.assertEqual(seed_id, "URS0000677DD8/1-71")

    # -----------------------------------------------------------------------

    def test_map_rnacentral_urs_wirh_db_accessions(self):

        db_accession = "hsa-mir-6753"

        urs_accession = rs.map_rnacentral_urs_wirh_db_accessions(db_accession)

        self.assertEqual(urs_accession, "URS000075CB3A_9606")

    # -----------------------------------------------------------------------

    def test_fetch_sequence_from_rnacentral(self):

        rnacentral_id = "URS000075CB3A"
        reference_seq = ''.join("CACCAGGGCAGAGCAGGGCUGAUCAUCUCACGUCA",
                                "GAGAGAGGGGAAGGGGCUGCCCAGUGAGCCCCCAC",
                                "AGGGCUCUACAUCUCCAGCUGGGCCUGGCUGGAGA",
                                "UCCCAGGGUCCCUGAAGGCCCCCGCCACCGUUCUG",
                                "GUCUGUCUCUGCCCUGGCACCCAG")

        sequence = rs.fetch_sequence_from_rnacentral(rnacentral_id)

        self.assertAlmostEqual(sequence, reference_seq)

    # -----------------------------------------------------------------------

    def test_map_sequence_segments(self):

        seed_sequence = "AUAGGUGCCGGUUUUUUAAAAAAACCCGACACCUUUAAUAUAUUAUAGGUGUCGGUAAAUUAAAAAUCGGCACCUAUAAUAAAUUAUAGGUGUCGGUUUUUUUAAAAACCGGCACCUGU"
        rnac_sequence = "AUAGGUGCCGGUUUUUUAAAAAAACCCGACACCUUUAAUAUAUUAUAGGUGUCGGUAAAUUAAAAAUCGGCACCUAUAAUAAAUUAUAGGUGUCGGUUUUUUUAAAAACCGGCACCUGU"
        new_seq = rs.map_sequence_segments(seed_sequence, rnac_sequence, no_segments=4)

        self.assertEqual(new_seq, "")

# -----------------------------------------------------------------------

if __name__ == '__main__':

    unittest.main()