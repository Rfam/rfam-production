import logging
import mysql.connector

from utils import RfamDB

rfam_search_url = "<https://rfam.org/family/{0}|{0}>"
pdb_search_url = "<https://www.rcsb.org/structure/{0}|{0}>"


def list_new_families():
    """
    List new families with 3D structures
    """
    conn = RfamDB.connect()
    cursor = conn.cursor()
    new_families_query = ("SELECT DISTINCT rfam_acc, pdb_id "
                          "FROM pdb_full_region "
                          "WHERE is_significant = 1 "
                          "AND rfam_acc NOT IN "
                          "(SELECT DISTINCT rfam_acc FROM pdb_full_region_old WHERE is_significant = 1);")
    try:
        cursor.execute("select count(distinct rfam_acc) from `pdb_full_region_old` where is_significant = 1;")
        print("Number of families with 3D before: {0}".format(cursor.fetchone()[0]))
        cursor.execute("select count(distinct rfam_acc) from `pdb_full_region` where is_significant = 1; ")
        print("Number of families with 3D after: {0}".format(cursor.fetchone()[0]))
        cursor.execute(new_families_query)
        new_families = cursor.fetchall()
        print("New families with 3D structures: {0}".format(new_families))
        for entry in new_families:
            rf_num = entry[0]
            pd_id = entry[1]
            print("Rfam accession number: {0}   PDB ID: {1}".
                  format(rfam_search_url.format(rf_num), pdb_search_url.format(pd_id)))
    except mysql.connector.Error as e:
        logging.debug("MySQL error has occurred: {0}".format(e))

    finally:
        cursor.close()
        conn.close()


if __name__ == '__main__':
    list_new_families()
