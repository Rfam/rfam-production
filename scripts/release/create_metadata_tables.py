import csv
import json
import os

import argparse

import mysql.connector

from utils import RfamDB
from config.gen_config import GENSEQ_FIELDS, RFAMSEQ_FIELDS, GENOME_FIELDS

genseq = []
rfamseq = []
genome = []


def generate_and_import(files_dir):
    """
    Make calls to generate the tables in the database, prepare and import the data
    :param files_dir: Directory location of the files with the metadata
    :return:
    """
    # for all the files in given directory
    for subdir, dirs, files in os.walk(files_dir):
        for f in files:
            filepath = os.path.join(subdir, f)
            if f.endswith(".jsonl"):
                # get the entries from the file
                genseq_entry, rfamseq_entry, genome_entry = get_entries_from_file(filepath)
                # add to main list
                genseq.append(genseq_entry)
                rfamseq.append(rfamseq_entry)
                genome.append(genome_entry)
    write_files(genseq, rfamseq, genome)
    create_tables()
    import_data_to_db()


def write_files(genseq, rfamseq, genome):
    """
    Write txt files ready for import to database
    :param genseq: list of entries
    :param rfamseq: list of entries
    :param genome: list of entries
    """
    with open('genseq_info_for_import.txt', 'w') as genseq_file:
        wr = csv.writer(genseq_file, delimiter='\t')
        wr.writerows(genseq)

    with open('rfamseq_info_for_import.txt', 'w') as rfamseq_file:
        wr = csv.writer(rfamseq_file, delimiter='\t')
        wr.writerows(rfamseq)

    with open('genome_info_for_import.txt', 'w') as genome_file:
        wr = csv.writer(genome_file, delimiter='\t')
        wr.writerows(genome)


def get_entries_from_file(json_file):
    """
    Get entries for the metadata tables from the .jsonl files
    :param json_file:
    :return: entry for each table
    """
    genseq_entry = []
    rfamseq_entry = []
    genome_entry = []
    with open(json_file, 'r') as f:
        data = json.load(f)
        upid = data["upid"]

        if data["genseq"]:
            genseq_entry.append(upid)
            for field in GENSEQ_FIELDS[1:]:
                if field in data["genseq"][0] and not data["genseq"][0][field] is None:
                    genseq_entry.append(data["genseq"][0][field])
                else:
                    genseq_entry.append('NULL')
        else:
            print('No genseq data for {0}, {1}'.format(json_file, upid))

        if data["rfamseq"]:
            for field in RFAMSEQ_FIELDS:
                if field in data["rfamseq"][0] and not data["rfamseq"][0][field] is None:
                    rfamseq_entry.append(data["rfamseq"][0][field])
                else:
                    rfamseq_entry.append('NULL')
        else:
            print('No rfamseq data for {0}, {1}'.format(json_file, upid))

        if data["genome"]:
            genome_entry.append(upid)
            for field in GENOME_FIELDS[1:]:
                if field in data["genome"] and not data["genome"][field] is None:
                    genome_entry.append(data["genome"][field])
                else:
                    genome_entry.append('NULL')
        else:
            print('No genome data for {0}, {1}'.format(json_file, upid))

    return genseq_entry, rfamseq_entry, genome_entry


def create_tables():
    """
    Create tables genome_temp, genseq_temp, rfamseq_temp ready to be populated with data
    """
    conn = RfamDB.connect()
    cursor = conn.cursor()
    try:
        cursor.execute("DROP TABLE IF EXISTS rfamseq_temp;")
        cursor.execute("CREATE TABLE rfamseq_temp LIKE rfamseq;")
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
    tables = ['genseq', 'rfamseq', 'genome']
    try:
        for table in tables:
            with open('{table}_info_for_import.txt'.format(table=table), 'r') as f:
                reader = csv.reader(f, delimiter='\t')
                for row in reader:
                    row = ['NULL' if val == '' else val for val in row]
                    row = [x.replace("'", "''") for x in row]
                    out = "'" + "', '".join(str(item) for item in row) + "'"
                    out = out.replace("'NULL'", 'NULL')
                    query = "INSERT INTO " + table + " VALUES (" + out + ")"
                    cursor.execute(query)
                conn.commit()

    except mysql.connector.Error as e:
        print("MySQL error has occurred: {0}".format(e))
        raise

    finally:
        cursor.close()
        conn.close()


def parse_args():
    """
    Parse the cli arguments
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-d', '--directory', help='location of files', required=True)
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    generate_and_import(args.directory)
