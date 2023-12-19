import argparse
import csv
import datetime
import logging
import os

import mysql.connector

from utils import RfamDB


def export_pdb_full_region_table(dest):
    """
    Export the data from pdb_full_region_table from rfam_live database into a text file
    """
    pdb_file = os.path.join(dest, "pdb_full_region_updated_" + str(datetime.date.today()) + ".txt")
    conn = RfamDB.connect()
    cursor = conn.cursor()
    try:
        cursor.execute("SELECT * FROM pdb_full_region ORDER BY rfam_acc;")
        results = cursor.fetchall()
        with open(pdb_file, 'w') as pdb_text:
            a = csv.writer(pdb_text, delimiter='\t')
            a.writerows(results)
    except mysql.connector.Error as e:
        logging.debug("MySQL error has occurred: {0}".format(e))
        raise

    finally:
        cursor.close()
        conn.close()


def parse_args():
    """
    Parse the cli arguments
    """
    parser = argparse.ArgumentParser(description='Create text file with pdb_full_region table contents')
    parser.add_argument('--dest-dir', help="destination directory to store output to", action="store")
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    export_pdb_full_region_table(dest=args.dest_dir)
