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

# Todo Need to generalize this to make the scripts independent and enable
# connecting to multiple databases simultaneously. Convert this to a class

# ---------------------------------IMPORTS-------------------------------------
import logging

import mysql.connector
from mysql.connector import errorcode

# Todo - NEED TO CLEAN THIS UP
# NEED TO CLEAN THIS UP AND PROFIDE A WAY
# from config.rfam_config import RFAMLIVEPUB  # rfam_live on public host
# from config.rfam_config import RFAMLIVE  # rfam_live on curation host
# from config.rfam_config import RFAMLIVELOC  # local instance of rfam_live

# from config.rfam_config import XFAMDEV
# from config.rfam_config import RFAMLOCAL
# from config.rfam_config import RFAMLIVE
# from config.rfam_config import RFAMLIVEPUB
# from config.rfam_config import RFAMREL

# RfamLive
from config.rfam_config import RFAMLIVE

# -----------------------------------------------------------------------------

# need to generalize this to enable DB setting upon implementation 
db_conf = RFAMLIVE

# -----------------------------------------------------------------------------


def connect(db_config=None):
    """
    Connects to a specific database and returns a mysql connection object.
    :param db_config: database config values
    :type db_config: dict or None
    :return: db connection
    """
    db_config = db_conf if db_config is None else db_config
    cnx = None
    try:
        cnx = mysql.connector.connect(user=db_config["user"],
                                      password=db_config["pwd"],
                                      host=db_config["host"],
                                      database=db_config["db"],
                                      port=db_config["port"])

    except mysql.connector.Error as err:

        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            print("Wrong username or password")
            logging.debug("Wrong username or password: {0}".format(err))

        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            print("Database does not exist")
            logging.debug("Database does not exist: {0}".format(err))

        else:
            print(err)
            logging.debug("MySQL error has occurred: {0}".format(err))

    return cnx


# -----------------------------------------------------------------------------


def disconnect(cnx):
    """
    Closes a database connection

    cnx: MySQL connection object
    """

    try:
        cnx.close()

    except Exception as e:
        print("Error closing database connection: {0}".format(e))
        logging.debug("Error closing database connection: {0}".format(e))

# -----------------------------------------------------------------------------
