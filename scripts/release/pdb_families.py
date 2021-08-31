import logging
import mysql.connector

from utils import RfamDB

rfam_search_url = "https://rfam.xfam.org/family/{0}"
pdb_search_url = "https://www.rcsb.org/structure/{0}"


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
        new_families_with_3D = cursor.fetchall()
        print("New families with 3D structures: {0}".format(new_families_with_3D))
        for entry in new_families_with_3D:
            print("RFAM accession number: {0}".format(entry[0]))
            print(rfam_search_url.format(entry[0]))
            print("PDB ID: {0}".format(entry[1]))
            print(pdb_search_url.format(entry[1]))
    except mysql.connector.Error as e:
        logging.debug("MySQL error has occurred: {0}".format(e))

    finally:
        cursor.close()
        conn.close()


if __name__ == '__main__':
    list_new_families()
