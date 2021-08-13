import logging
import mysql.connector
from config.rfam_config import RFAMLIVE
from utils import RfamDB

db_conf = RFAMLIVE


def list_new_families():
    """
    List new families with 3D structures
    """

    # when ready to apply changes the pdb_full_region will be the 'new'

    conn = RfamDB.connect()
    cursor = conn.cursor()
    # new_families_query = ("SELECT DISTINCT rfam_acc "
    #                       "FROM pdb_full_region "
    #                       "WHERE is_significant = 1 "
    #                       "AND rfam_acc NOT IN "
    #                       "(SELECT DISTINCT rfam_acc FROM pdb_full_region_temp WHERE is_significant = 1);")
    new_families_query = ("SELECT DISTINCT rfam_acc "
                          "FROM pdb_full_region_temp "
                          "WHERE is_significant = 1 "
                          "AND rfam_acc NOT IN "
                          "(SELECT DISTINCT rfam_acc FROM pdb_full_region WHERE is_significant = 1);")
    try:
        # cursor.execute('select count(distinct rfam_acc) from `pdb_full_region_old` where is_significant = 1;')
        cursor.execute('select count(distinct rfam_acc) from pdb_full_region where is_significant = 1;')
        print('Number of families with 3D before: {0}'.format(cursor.fetchall()))
        # cursor.execute('select count(distinct rfam_acc) from `pdb_full_region` where is_significant = 1; ')
        cursor.execute('select count(distinct rfam_acc) from pdb_full_region_temp where is_significant = 1; ')
        print('Number of families with 3D after: {0}'.format(cursor.fetchall()))
        cursor.execute(new_families_query)
        print('New families with 3D structures: {0}'.format(cursor.fetchall()))
        conn.commit()
    except mysql.connector.Error as e:
        logging.debug("MySQL error has occurred: {0}".format(e))

    finally:
        cursor.close()
        conn.close()


if __name__ == '__main__':
    list_new_families()
