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

'''
TO DO: Need to generalize this to make the scripts independent and enable
       connecting to multiple databases simultaneously. Convert this to a class
'''

# ---------------------------------IMPORTS-------------------------------------
import mysql.connector
from mysql.connector import errorcode

from config.rfam_config import RFAMLIVEPUB  # rfam_live on public host
#from config.rfam_config import RFAMLIVE  # rfam_live on curation host
from config.rfam_config import RFAMLIVELOC  # local instance of rfam_live

from config.rfam_config import XFAMDEV
from config.rfam_config import RFAMLOCAL
from config.rfam_config import RFAMLIVE
from config.rfam_config import RFAMLIVEPUB
from config.rfam_config import RFAMREL
# -----------------------------------------------------------------------------

# need to generalize this to enable DB setting upon implementation 
db_conf = RFAMLIVE

# -----------------------------------------------------------------------------


def connect():
    """
    Connects to a specific database and returns a mysql connection object
    """

    cnx=None
    try:
        cnx = mysql.connector.connect(user=db_conf["user"],
                                      password=db_conf["pwd"],
                                      host=db_conf["host"],
                                      database=db_conf["db"],
                                      port=db_conf["port"])

    except mysql.connector.Error as err:

        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            print "Wrong username or password"

        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            print "Database does not exist"

        else:
            print err

    return cnx

# -----------------------------------------------------------------------------


def disconnect(cnx):
    """
    Closes a database connection

    cnx: MySQL connection object
    """

    try:
        cnx.close()

    except:
        print "Error closing database connection"

# -----------------------------------------------------------------------------
