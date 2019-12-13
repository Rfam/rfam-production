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

# -----------------------------------------------------------------------

if __name__ == '__main__':
	
	unittest.main()
	
