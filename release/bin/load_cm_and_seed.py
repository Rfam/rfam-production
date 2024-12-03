"""
Update CM and seed alignment blobs in the Rfam MySQL database.

The blobs are used by the following endpoints:
https://preview.rfam.org/family/RFXXXXX/cm
https://preview.rfam.org/family/RFXXXXX/alignment?acc=RFXXXXX&format=stockholm&download=0

Usage:
python load_cm_seed_in_db.py /path/to/Rfam.seed /path/to/Rfam.cm
"""

import os
import sys
import tempfile

import pymysql.cursors

# from utils import RfamDB
# from utils.getters import get_all_rfam_accessions


def get_all_rfam_accessions():
    """
    Fetch a list of all Rfam families from the SVN repository.
    """
    rfam_accessions = []
    svn_url = "https://svn.rfam.org/svn/data_repos/trunk/Families/"
    svn_list = tempfile.NamedTemporaryFile()
    cmd = "svn list {} > {}".format(svn_url, svn_list.name)
    os.system(cmd)
    with open(svn_list.name, "r") as f_svn_list:
        for line in f_svn_list:
            if line.startswith("RF"):
                rfam_accessions.append(line.strip().replace("/", ""))
    print("Found {} accessions on SVN".format(len(rfam_accessions)))
    return rfam_accessions


def load_cm_seed_in_db(rfam_acc, seed_file, cm_file, cursor, cnx):
    seed_gz_file = "{0}.seed.gz".format(rfam_acc)
    cmd = "esl-afetch {0} {1} | gzip > {2}".format(seed_file, rfam_acc, seed_gz_file)
    os.system(cmd)
    with open(seed_gz_file, "rb") as f:
        seed_gzip = f.read()
    os.remove(seed_gz_file)

    cm_gz_file = "{}.cm.gz".format(rfam_acc)
    cmd = "cmfetch {0} {1} | gzip > {2}".format(cm_file, rfam_acc, cm_gz_file)
    os.system(cmd)
    with open("{0}.cm.gz".format(rfam_acc), "rb") as f:
        cm_gzip = f.read()
    os.remove(cm_gz_file)

    query = "REPLACE INTO _annotated_file (rfam_acc, seed, cm) VALUES (%s, %s, %s)"
    cursor.execute(query, (rfam_acc, seed_gzip, cm_gzip))
    cnx.commit()


def main(seed_file, cm_file):
    db_config = {
        "user": "admin",
        "pwd": "kQ02PjvA",
        "host": "mysql-rfam-rel.ebi.ac.uk",
        "db": "rfam_live",
        "port": 4442,
    }
    cnx = pymysql.connect(
        user=db_config["user"],
        password=db_config["pwd"],
        host=db_config["host"],
        database=db_config["db"],
        port=db_config["port"],
    )
    cursor = cnx.cursor()
    for rfam_acc in get_all_rfam_accessions():
        print(rfam_acc)
        load_cm_seed_in_db(rfam_acc, seed_file, cm_file, cursor, cnx)
    cursor.close()
    print("Done")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise Exception("Specify seed and CM file arguments")

    if not sys.argv[1] or not os.path.exists(sys.argv[1]):
        raise Exception("Error: seed file {} does not exist".format(sys.argv[1]))

    if not sys.argv[2] or not os.path.exists(sys.argv[2]):
        raise Exception("Error: CM file {} does not exist".format(sys.argv[2]))

    main(sys.argv[1], sys.argv[2])
