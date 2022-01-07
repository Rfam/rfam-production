#!/usr/bin/env python

from contextlib import contextmanager
import tempfile
from ftplib import FTP
import subprocess as sp
from datetime import date

from config import rfam_config

QUERY = """
select
    accession, 
    rfam_acc, 
    rfam_id
from full_region
join family using (rfam_acc)
join rfamseq using (rfamseq_acc)
"""


@contextmanager
def query():
    db = rfam_config.RFAMLIVE
    with tempfile.NamedTemporaryFile() as tmp:
        process = sp.Popen(
            args=[
                "mysql",
                "--host",
                db["host"],
                "--port",
                str(db["port"]),
                "--user",
                db["user"],
                "--password=%s" % db["pwd"],
                db["db"],
            ],
            stdin=sp.PIPE,
            stdout=tmp,
        )

        (stdout, stderr) = process.communicate(input=QUERY.encode('utf-8'))
        assert stdout is None
        assert stderr is None
        yield tmp.name


def upload_file(filename):
    config = rfam_config.ENA_FTP
    with FTP(config['hostname'], config['user'], config['password']) as ftp, open(
        filename, "rb"
    ) as raw:
        now = date.today().isoformat()
        filename = 'Rfam_%s_rfam2embl_crossrefs.txt' % now
        ftp.storbinary('STOR %s' % filename, raw)


def main():
    with query() as temp_file:
        upload_file(temp_file)


if __name__ == "__main__":
    main()
