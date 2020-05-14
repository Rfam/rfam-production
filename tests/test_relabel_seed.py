import unittest
import scripts.preprocessing.relabel_seed as rs

# -----------------------------------------------------------------------

class TestRelabelSeedFunctions(unittest.TestCase):

    def test_load_fasta_file_to_dict(self):
        pass


    def test_sequence_to_md5(self):

        "https://rnacentral.org/api/v1/rna/URS0000677DD8.json"

        # RNAcentral accession: URS0000677DD8
        sequence = "UCAAUAAUGAAAUCUUCUGAUUUGGUGAGAAAUAAUGCCUUAAAAUUACACUCAAUAGGAUUAUGCUGAGG"

        md5_hash = rs.sequence_to_md5(sequence)

        self.assertEqual(md5_hash, "e25b1f5967d21e39e394edb4d3054d0e")


    def test_fetch_RNAcentral_id(self):
        sequence = "UCAAUAAUGAAAUCUUCUGAUUUGGUGAGAAAUAAUGCCUUAAAAUUACACUCAAUAGGAUUAUGCUGAGG"

        rnacentral_id = rs.fetch_RNAcentral_id(sequence)

        self.assertEqual(rnacentral_id, "URS0000677DD8")

# -----------------------------------------------------------------------

if __name__ == '__main__':

    unittest.main()