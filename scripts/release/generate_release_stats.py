import csv
import logging
import sys

import mysql.connector

from utils import RfamDB


def update_file(release_stats_file):
    """
    Update the release_stats.tsv file in the FTP folder to include the latest version number, release date, and
    number of families
    """
    conn = RfamDB.connect()
    cursor = conn.cursor()
    try:
        version = cursor.execute("select rfam_release from `version`;")
        release_date = cursor.execute("select rfam_release_date from `version`;")
        num_families = cursor.execute("select count(distinct rfam_acc) from `family`;")
    except mysql.connector.Error as e:
        logging.debug("MySQL error has occurred: {0}".format(e))
        raise e
    except Exception as e:
        raise e
    finally:
        cursor.close()
        conn.close()

    with open(release_stats_file, 'w') as rsf:
        writer = csv.writer(rsf, delimiter='\t')
        row = [version, release_date, num_families]
        writer.writerow(row)


if __name__ == '__main__':
    stats_file = sys.argv[1:]
    update_file(stats_file)
