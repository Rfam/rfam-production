import argparse
import logging
import mysql.connector

from config.rfam_config import RFAMREL
from utils import RfamDB


def update_table(pdb_text_file):
    """
    Update pdb_full_region table in RfamRel database.
    """
    conn = RfamDB.connect(db_config=RFAMREL)
    cursor = conn.cursor()
    try:
        cursor.execute("CREATE TABLE pdb_full_region_temp LIKE pdb_full_region;")
        cursor.execute("LOAD DATA LOCAL INFILE '{0}' INTO TABLE pdb_full_region_temp "
                       "COLUMNS TERMINATED BY '\t';".format(pdb_text_file))
        cursor.execute("DROP TABLE pdb_full_region_old;")
        cursor.execute("RENAME TABLE pdb_full_region TO pdb_full_region_old,pdb_full_region_temp TO pdb_full_region;")
    except mysql.connector.Error as e:
        logging.debug("MySQL error has occurred: {0}".format(e))
    finally:
        cursor.close()
        conn.close()


def parse_args():
    """
    Parse the cli arguments when calling this script to insert a text file to the PDB table in the database
    """
    parser = argparse.ArgumentParser(description='Create PDB full region table in RfamRel and import new data')
    parser.add_argument('-f', '--file', help='Text file with data to import to pdb_full_region_temp', required=True)
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    if args.file:
        update_table(args.file)
    else:
        print("Please provide a text file to import data to pdb_full_region_temp")
