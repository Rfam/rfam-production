'''
Created on 18 Feb 2016

@author: ikalvari

TO DO: Need to generalize this to make the scripts independent and enable
       connecting to multiple databases simultaneously. Convert this to a class
'''

# ---------------------------------IMPORTS-------------------------------------
import mysql.connector
from mysql.connector import errorcode
from config.rfam_config import RFAMLIVEPUB  # rfam_live on public host
from config.rfam_config import RFAMLIVE  # rfam_live on curation host
from config.rfam_config import RFAMLIVELOC  # local instance of rfam_live

# -----------------------------------------------------------------------------
# need to generalize this enabling database selection upon execution
db_conf = RFAMLIVE

# -----------------------------------------------------------------------------


def connect():
    '''
    Connects to a specific database and returns a mysql connection object.

    '''

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
    '''
        Closes a database connection.

        cnx: MySQL connection object
    '''

    try:
        cnx.close()

    except:
        print "Error closing database connection"

# -----------------------------------------------------------------------------
