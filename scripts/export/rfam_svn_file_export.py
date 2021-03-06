"""
Copyright [2009-2017] EMBL-European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
     http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""
"""
Description: Export script fetching family specific files from the SVN 
             (e.g. "SEED", "CM")
"""
# ---------------------------------IMPORTS-------------------------------------

import os
import sys
import subprocess
import shutil
import string
from utils import RfamDB

# -----------------------------------------------------------------------------

FILE_TYPES = ["CM", "SEED"]
SVN_CHECKOUT = "rfco %s"


# -----------------------------------------------------------------------------


def export_rfam_family_files(f_types, out_dir):
    """
    Fetches all Rfam family accessions from rfam_live, checks out each
    family and copies the files in f_types in their corresponding
    directories

    f_types: A list of file type keywords we need to
             export (e.g. ["SEED", "CM"])
    out_dir: The path to the output directory. If it does not exist it will
             be created
    """

    # Create the output directory if it does not exist
    if (not os.path.exists(out_dir)):
        os.mkdir(out_dir)

    # if current working directory isn't out_dir, change directory
    if (string.find(os.getcwd(), out_dir) == -1):
        os.chdir(out_dir)

    # generate specific output directories for each file type
    file_path = ''
    for f_type in f_types:
        file_path = os.path.join(out_dir, f_type)
        if (not os.path.exists(file_path)):
            os.mkdir(file_path)
        file_path = ''

    # get DB connection handle
    cnx = RfamDB.connect()

    # get mysql cursor
    cursor = cnx.cursor(buffered=True)

    # execute query
    cursor.execute("SELECT rfam_acc FROM family")

    cmd = ''

    # fetch files for all Rfam family accessions
    for rfam_acc in cursor:

        rfam_acc = str(rfam_acc[0])
        cmd = SVN_CHECKOUT % rfam_acc

        # Check out family in out_dir using rfco on lsf
        subprocess.call(cmd, shell=True)

        # path to
        fam_dir = os.path.join(out_dir, rfam_acc)

        # copy files and rename
        for f_type in f_types:
            filename = rfam_acc + '.' + f_type.lower()
            """
            if (f_type == "SEED"):
                # 1. open out file handler
                seed_out_fp = open(
                    os.path.join(os.path.join(out_dir, f_type), filename), 'w')
                # 2. open desc handler
                desc_fp = open(os.path.join(fam_dir, "DESC"), 'r')
                # 3. write desc to outfile
                seed_out_fp.writelines(desc_fp.readlines())
                seed_out_fp.write('\n')
                desc_fp.close()
                # 4. open seed and write in outfile

                continue
            """
            shutil.copyfile(os.path.join(fam_dir, f_type),
                            os.path.join(os.path.join(out_dir, f_type), filename))

        # delete family dir
        shutil.rmtree(fam_dir)

        filename = ''
        fam_dir = ''
        cmd = ''

    # close DB connection
    cursor.close()
    RfamDB.disconnect(cnx)


# -----------------------------------------------------------------------------


def usage():
    """
    Displays information on how to run rfam_svn_file_export
    """

    print "\nUsage:\n------"

    print "\npython rfam_svn_file_export.py out_dir"
    print "\nout_dir: Path to an output directory"
    print "\n-h option for usage\n"


# -----------------------------------------------------------------------------

if __name__ == '__main__':

    out_dir = sys.argv[1]

    # minor input checks
    if os.path.isdir(out_dir):
        export_rfam_family_files(FILE_TYPES, out_dir)

    else:
        usage()
