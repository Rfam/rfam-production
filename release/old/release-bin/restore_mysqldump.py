import argparse
import logging

import mysql.connector
from config.rfam_config import RFAMLIVEPUB, RFAMREL
from utils import RfamDB

DB_CONFIG = None


def restore_dump():
    """
    Connect to the given database, create a new schema, and restore the database
    """
    conn = RfamDB.connect(db_config=DB_CONFIG)
    cursor = conn.cursor()
    try:
        cursor.execute("Create schema rfam_X_Y;")
        cursor.execute("Use rfam_X_Y;")
        cursor.execute("source rfam_live_relX.sql")
    except mysql.connector.Error as e:
        print("MySQL error has occurred: {0}".format(e))
        raise
    finally:
        cursor.close()
        conn.close()


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rel", help="Restore mysqldump for RfamRel")
    parser.add_argument("--public", help="Restore mysqldump for Rfam Public")

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()
    if args.rel:
        DB_CONFIG = RFAMREL
    if args.public:
        DB_CONFIG = RFAMLIVEPUB
    restore_dump()
