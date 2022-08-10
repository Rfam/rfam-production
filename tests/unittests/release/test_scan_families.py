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

    @patch('scan_families.subprocess.call')
    def test_checkout_family(self, mock_subprocess):
        scan_families.checkout_family('RF00001')
        mock_subprocess.assert_called_with("rfco.pl RF00001", shell=True)

    @patch('scan_families.subprocess.call')
    def test_submit_new_rfsearch_job_rfmake_true(self, mock_subprocess):
        rfmake_true_cmd = 'bsub -M 2000 -R "rusage[mem=2000]" -o /family/dir/path/auto_rfsearch.out -e ' \
                          '/family/dir/path/auto_rfsearch.err -q bigmem "cd /family/dir/path && ' \
                          'rfsearch.pl -scpu 0 -cnompi -q short && rfmake.pl -local"'
        scan_families.submit_new_rfsearch_job('/family/dir/path', rfmake=True)
        mock_subprocess.assert_called_with(rfmake_true_cmd, shell=True)

    @patch('scan_families.subprocess.call')
    def test_submit_new_rfsearch_job_rfmake_false(self, mock_subprocess):
        rfmake_false_cmd = 'bsub -M 2000 -R "rusage[mem=2000]" -o /family/dir/path/auto_rfsearch.out -e ' \
                          '/family/dir/path/auto_rfsearch.err -q bigmem "cd /family/dir/path && ' \
                          'rfsearch.pl -scpu 0 -cnompi -q short"'
        scan_families.submit_new_rfsearch_job('/family/dir/path')
        mock_subprocess.assert_called_with(rfmake_false_cmd, shell=True)
