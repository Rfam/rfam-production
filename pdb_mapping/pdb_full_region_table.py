import argparse
import csv
import logging

import mysql.connector

from pdb_mapping.exceptions import CheckRowsError, NoFileGivenError
from utils import RfamDB
from config.rfam_config import RFAMREL

DB_CONFIG = None


def create_pdb_temp_table(pdb_file):
    """
    Create the pdb_full_region_temp table and populate with data from the pdb text file.
    :param pdb_file: Text file with data to import to pdb_full_region_temp
    """
    conn = RfamDB.connect(db_config=DB_CONFIG)
    cursor = conn.cursor()
    try:
        cursor.execute("DROP TABLE IF EXISTS pdb_full_region_temp;")
        cursor.execute("CREATE TABLE pdb_full_region_temp LIKE pdb_full_region;")
        with open(pdb_file) as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                cursor.execute("""INSERT INTO pdb_full_region_temp(rfam_acc, pdb_id, chain, pdb_start,
                                pdb_end, bit_score, evalue_score, cm_start, cm_end, hex_colour, is_significant)
                                VALUES(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s) """, row)
                conn.commit()
    except mysql.connector.Error as e:
        logging.debug("MySQL error has occurred: {0}".format(e))
        raise

    finally:
        cursor.close()
        conn.close()


def qc_checks():
    """
    Execute quality control checks before we update the table
    """
    conn = RfamDB.connect(db_config=DB_CONFIG)
    cursor = conn.cursor()
    try:
        cursor.execute("SELECT COUNT(*) FROM pdb_full_region_temp;")
        num_rows_pdb_temp = cursor.fetchone()[0]
        cursor.execute("SELECT COUNT(*) FROM pdb_full_region;")
        num_rows_pdb = cursor.fetchone()[0]
        try:
            if num_rows_pdb_temp - num_rows_pdb > 100:
                raise CheckRowsError(num_rows_pdb, num_rows_pdb_temp)
            if num_rows_pdb_temp < 10000:
                raise CheckRowsError(num_rows_pdb, num_rows_pdb_temp)
        except CheckRowsError:
            logging.exception("Do not rename tables until error is resolved.")
            raise
        except Exception as e:
            logging.exception(e)

    except mysql.connector.Error as e:
        logging.debug("MySQL error has occurred: {0}".format(e))
        raise 

    finally:
        cursor.close()
        conn.close()


def rename_pdb_table():
    """
    Rename pdb_full_region to pdb_full_region_old and rename pdb_full_region_temp to pdb_full_region.
    """
    conn = RfamDB.connect(db_config=DB_CONFIG)
    cursor = conn.cursor()
    try:
        cursor.execute("DROP TABLE IF EXISTS pdb_full_region_old;")
        cursor.execute("RENAME TABLE pdb_full_region TO pdb_full_region_old, "
                       "pdb_full_region_temp TO pdb_full_region;")
    except mysql.connector.Error as e:
        logging.debug("MySQL error has occurred: {0}".format(e))
        raise
    finally:
        cursor.close()
        conn.close()


def parse_args():
    """
    Parse the cli arguments when calling this script to insert a text file to the PDB table in the database.
    """
    parser = argparse.ArgumentParser(description='Create PDB full region table and import new data')
    parser.add_argument('-f', '--file', help='Text file with data to import to pdb_full_region_temp', required=True)
    parser.add_argument('-db', '--database', help='Specify which database config values to use ', required=False)
    parser.add_argument('-nqc', '--no_quality_control_checks',
                        help='Whether to execute QC checks ', action='store_true', required=False)
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    if args.database == 'rfam-rel':
        DB_CONFIG = RFAMREL
    if args.file and args.no_quality_control_checks:
        create_pdb_temp_table(args.file)
        rename_pdb_table()
    elif args.file:
        create_pdb_temp_table(args.file)
        qc_checks()
        rename_pdb_table()
    else:
        raise NoFileGivenError("Please provide a text file to import data to pdb_full_region_temp")
