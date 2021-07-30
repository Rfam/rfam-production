import logging
import mysql.connector
import csv
from config.rfam_config import RFAMLIVE

db_conf = RFAMLIVE


# merge with RfamDB?
def connect():
    """
    Connect to the database and return the connection
    """
    try:
        conn = mysql.connector.connect(user=db_conf["user"],
                                       password=db_conf["pwd"],
                                       host=db_conf["host"],
                                       database=db_conf["db"],
                                       port=db_conf["port"])
        return conn
    except mysql.connector.Error as e:
        logging.debug("MySQL error has occurred: {0}".format(e))


def create_pdb_temp_table():
    """
    Create the pdb_full_region_temp table and populate with data from the pdb text file
    """

    conn = connect()
    cursor = conn.cursor()
    try:
        cursor.execute('CREATE TABLE IF NOT EXISTS pdb_full_region_temp LIKE pdb_full_region')
        with open('pdb_full_region_2021-07-20.txt') as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                cursor.execute("""INSERT INTO pdb_full_region_temp(rfam_acc, pdb_id, chain, pdb_start,
                pdb_end, bit_score, evalue_score, cm_start, cm_end, hex_colour, is_significant)
                VALUES(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s) """, row)
    except mysql.connector.Error as e:
        logging.debug("MySQL error has occurred: {0}".format(e))

    finally:
        cursor.close()
        conn.close()
        # conn.commit() will automatically be called when Python leaves the outer `with` statement
        # Neither crs.close() nor conn.close() will be called upon leaving the `with` statement!!
        # TODO context mgr


def rename_pdb_table():
    """
    Rename pdb_full_region to pdb_full_region_old and rename pdb_full_region_temp to pdb_full_region
    """
    conn = connect()
    cursor = conn.cursor()
    try:
        # RENAME TABLE, unlike ALTER TABLE, can rename multiple tables within a single statement
        cursor.execute('RENAME TABLE pdb_full_region TO pdb_full_region_old, '
                       'pdb_full_region_temp TO pdb_full_region;')
        print('complete')
    except mysql.connector.Error as e:
        logging.debug("MySQL error has occurred: {0}".format(e))
    finally:
        cursor.close()
        conn.close()


if __name__ == '__main__':
    create_pdb_temp_table()
    #  don't call until sure of change - will overwrite pdb_full_region
    # rename_pdb_table()
