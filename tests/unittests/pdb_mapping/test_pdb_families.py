import unittest
from mock import patch

from pdb_mapping.pdb_families import list_new_families


class PdbTest(unittest.TestCase):

    @patch('pdb_mapping.pdb_families.RfamDB.connect')
    def test_list_new_families_is_successful(self, mock_connect):
        new_families_query = ("SELECT DISTINCT rfam_acc, pdb_id "
                              "FROM pdb_full_region "
                              "WHERE is_significant = 1 "
                              "AND rfam_acc NOT IN "
                              "(SELECT DISTINCT rfam_acc FROM pdb_full_region_old WHERE is_significant = 1);")
        conn = mock_connect.return_value
        cursor = conn.cursor.return_value
        list_new_families()
        self.assertEqual(5, cursor.execute.call_count)
        cursor.execute.assert_called_with(new_families_query)

    @patch('pdb_mapping.pdb_families.RfamDB.connect')
    def test_list_new_families_raises_mysql_error(self, mock_connect):
        conn = mock_connect.return_value
        cursor = conn.cursor.return_value
        cursor.execute.side_effect = Exception
        with self.assertRaises(Exception):
            list_new_families()

    @patch('pdb_mapping.pdb_families.RfamDB.connect')
    def test_list_new_families_raises_exception(self, mock_connect):
        conn = mock_connect.return_value
        cursor = conn.cursor.return_value
        cursor.execute.side_effect = Exception
        with self.assertRaises(Exception):
            list_new_families()


if __name__ == '__main__':
    unittest.main()
