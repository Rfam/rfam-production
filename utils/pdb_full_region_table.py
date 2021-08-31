import argparse
import csv
import logging

import mysql.connector

from utils import RfamDB


def create_pdb_temp_table(pdb_file):
    """
    Create the pdb_full_region_temp table and populate with data from the pdb text file
    """

    conn = RfamDB.connect()
    cursor = conn.cursor()
    try:
        cursor.execute("CREATE TABLE pdb_full_region_temp LIKE pdb_full_region")
        with open(pdb_file) as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                cursor.execute("""INSERT INTO pdb_full_region_temp(rfam_acc, pdb_id, chain, pdb_start,
                                pdb_end, bit_score, evalue_score, cm_start, cm_end, hex_colour, is_significant)
                                VALUES(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s) """, row)
                conn.commit()
    except mysql.connector.Error as e:
        logging.debug("MySQL error has occurred: {0}".format(e))

    finally:
        cursor.close()
        conn.close()


def rename_pdb_table():
    """
    Rename pdb_full_region to pdb_full_region_old and rename pdb_full_region_temp to pdb_full_region
    """
    conn = RfamDB.connect()
    cursor = conn.cursor()
    try:
        cursor.execute("DROP TABLE pdb_full_region_old;")
        cursor.execute("RENAME TABLE pdb_full_region TO pdb_full_region_old, "
                       "pdb_full_region_temp TO pdb_full_region;")
    except mysql.connector.Error as e:
        logging.debug("MySQL error has occurred: {0}".format(e))
    finally:
        cursor.close()
        conn.close()


def parse_args():
    """
    Parse the cli arguments when calling this script to insert a text file to the PDB table in the database
    """
    parser = argparse.ArgumentParser(description='Create PDB full region table and import new data')
    parser.add_argument('-f', '--file', help='Text file with data to import to pdb_full_region_temp', required=True)
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    if args.file:
        create_pdb_temp_table(args.file)
        rename_pdb_table()
    else:
        print("Please provide a text file to import data to pdb_full_region_temp")
