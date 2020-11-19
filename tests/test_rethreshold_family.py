import unittest

import scripts.release.rethreshold_family as rf

# -----------------------------------------------------------------------

class TestFamilyRethresholdingFunctions(unittest.TestCase):
    def test_extract_unique_seeds_from_seedoutlist(self):
        """
        """

        seedoutlist = "../data/RF02058_seedoutlist.txt"

        seed_accs = rf.extract_unique_seeds_from_seedoutlist(seedoutlist)

        # The function should be able to extract the 2 sequences in the
        # seedoutlist file
        self.assertEqual(len(seed_accs.keys()), 2)

    # ---------------------------------------------------------------------

    def test_checkout_family(self):
        """
        Tests if rfco.pl command works and the Rfam family is checked out
        successfully
        """
        rfam_acc = 'RF02567'

        rf.checkout_family(rfam_acc)

        self.assertExists()

    # ---------------------------------------------------------------------

    def test_extract_scores_from_outlist_file(self):
        """

        :return:
        """

        outlist = "../data/outlist"
        scores = rf.extract_scores_from_outlist_file(outlist)

        self.assertDictEqual()

    # ---------------------------------------------------------------------

if __name__ == '__main__':
    unittest.main()
