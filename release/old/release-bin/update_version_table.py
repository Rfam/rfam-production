import argparse
import datetime
import logging

import mysql.connector
from utils import RfamDB


def update_version(release_version):
    """
    Update the fields in the version table for the latest release
    """
    conn = RfamDB.connect()
    cursor = conn.cursor()
    num_families_query = "select count(*) from family;"
    version_query = (
        "UPDATE version SET rfam_release = {release}, rfam_release_date = '{date}', number_families = {"
        "num_families}, embl_release = 132;"
    )
    try:
        today_date = str(datetime.date.today())
        cursor.execute(num_families_query)
        num_families = cursor.fetchone()[0]
        cursor.execute(
            version_query.format(
                release=release_version, date=today_date, num_families=num_families
            )
        )
        conn.commit()
    except mysql.connector.Error as e:
        logging.debug("MySQL error has occurred: {0}".format(e))
        raise e
    except Exception as e:
        raise e

    finally:
        cursor.close()
        conn.close()


def parse_args():
    """
    Basic argument parsing using Python's argparse

    return: Argparse parser object
    """

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--version", help="Release version number, e.g. 14.8", required=True
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    update_version(args.version)
