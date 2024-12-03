import csv
import datetime

import mysql.connector
from utils import RfamDB

genseq = []
rfamseq = []
genome = []


def generate_and_import():
    """
    Make calls to generate the tables in the database and import the data
    :return:
    """
    create_tables()
    import_data_to_db()


def create_tables():
    """
    Create tables genome_temp, genseq_temp, rfamseq_temp ready to be populated with data
    """
    conn = RfamDB.connect()
    cursor = conn.cursor()
    try:
        cursor.execute("DROP TABLE IF EXISTS rfamseq_temp;")
        cursor.execute("CREATE TABLE rfamseq_temp LIKE rfamseq;")
        cursor.execute("ALTER TABLE rfamseq_temp MODIFY description LONGTEXT;")
        cursor.execute("DROP TABLE IF EXISTS genseq_temp;")
        cursor.execute("CREATE TABLE genseq_temp LIKE genseq;")
        cursor.execute("DROP TABLE IF EXISTS genome_temp;")
        cursor.execute("CREATE TABLE genome_temp LIKE genome;")
    except mysql.connector.Error as e:
        print("MySQL error has occurred: {0}".format(e))
        raise
    finally:
        cursor.close()
        conn.close()


def import_data_to_db():
    """
    Import data from the text files to the rfamseq, genseq, genome tables
    """
    conn = RfamDB.connect()
    cursor = conn.cursor()
    tables = ["genseq", "rfamseq", "genome"]
    try:
        for table in tables:
            with open("{table}_info_for_import.txt".format(table=table), "r") as f:
                reader = csv.reader(f, delimiter="\t")
                for row in reader:
                    if len(row) == 0:
                        print("Empty row in ", table)
                    else:
                        row = ["NULL" if val == "" else val for val in row]
                        row = [x.replace("'", "''") for x in row]
                        if table == "genome":
                            # set is_reference, is_representative to 0
                            row[18] = row[19] = "0"
                            # set created, updated to today's date
                            entry_date = datetime.datetime.now().strftime(
                                "%Y-%m-%d %H:%M:%S"
                            )
                            row[20] = row[21] = entry_date
                        out = "'" + "', '".join(str(item) for item in row) + "'"
                        out = out.replace("'NULL'", "NULL")
                        query = "INSERT INTO " + table + "_temp VALUES (" + out + ")"
                        cursor.execute(query)
                conn.commit()

    except mysql.connector.Error as e:
        print("MySQL error has occurred: {0}".format(e))
        raise

    finally:
        cursor.close()
        conn.close()


if __name__ == "__main__":
    generate_and_import()
