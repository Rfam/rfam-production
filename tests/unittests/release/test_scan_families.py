import unittest
from mock import patch, Mock

import scan_families


class TestScans(unittest.TestCase):

    def setUp(self):
        self.test_file = "data/rfam_accs.txt"
        self.accs_list = ["RF00001", "RF00002", "RF00005", "RF04192"]

    @patch('scan_families.checkout_and_search_family')
    @patch('scan_families.list_rfam_accessions_from_file')
    def test_run_searches(self, mock_list_rfam_accs, mock_checkout):
        mock_list_rfam_accs.return_value = self.accs_list
        scan_families.run_searches(self.test_file, Mock())
        self.assertEqual(4, mock_checkout.call_count)

    def test_load_rfam_accessions_from_file(self):
        self.assertEqual(scan_families.list_rfam_accessions_from_file(self.test_file), self.accs_list)

    @patch('scan_families.os')
    @patch('scan_families.checkout_family')
    @patch('scan_families.submit_new_rfsearch_job')
    def test_checkout_and_search_family(self, mock_submit, mock_checkout, mock_os, ):
        mock_os.path.exists.side_effect = False, True
        scan_families.checkout_and_search_family('RF00001', Mock())
        self.assertEqual(1, mock_submit.call_count)
        self.assertEqual(1, mock_checkout.call_count)
