'''
Created on 20 Apr 2016

@author: ikalvari
'''

# ---------------------------------IMPORTS-------------------------------------

import os
import sys
from utils import RfamDB

# -----------------------------------------------------------------------------


def check_ss_images(cnx_obj, no_fams):
    '''
        Checks that secondary structure images have been generated for all
        new families and returns True or False accordingly.
    '''

    query = ("SELECT ss.type,count(*)\n"
             "FROM secondary_structure_image ss, family f\n"
             "WHERE f.rfam_acc=ss.rfam_acc\n"
             "AND ss.image is not NULL\n"
             "AND f.author like \'%Arga%\'\n"
             "GROUP BY ss.type")

    # TO BE IMPLEMENTED

# -----------------------------------------------------------------------------


def check_sunburst(cnx_obj):
    '''
        Looks up sunburst table and checks that there're entries for all
        families in the view process.
    '''
    # families number or get that dynamically

    query = ("SELECT s.type,count(*)\n"
             "FROM sunburst s, family f\n"
             "WHERE f.rfam_acc=s.rfam_acc\n"
             "AND s.data is not NULL\n"
             "AND f.author like \'%Arga%\'\n"
             "GROUP BY s.type")

    # TO BE IMPLEMENTED

# -----------------------------------------------------------------------------


def count_rchie_diagrams(cnx_obj, no_fams):
    '''
        Counts the number of rchie diagrams generated
    '''

    query = ("SELECT count(*) from secondary_structure_image\n"
             "WHERE type=\'rchie\' and image is not NULL\n")

    # TO BE IMPLEMENTED


# -----------------------------------------------------------------------------


def check_alignment_and_tree(cnx_obj, no_fams):
    '''
        Checks all types of files for the new families in alignment_and_tree
        table.
    '''

    query = ("SELECT ant.type, count(*)\n"
             "FROM alignment_and_tree ant, family f\n"
             "WHERE f.rfam_acc=ant.rfam_acc\n"
             "AND f.author like \'%Arga%\'\n"
             "GROUP BY ant.type")

    # TO BE IMPLEMENTED

# -----------------------------------------------------------------------------


def check_html_alignment(cnx_obj, no_fams):
    '''
        Checks if there're entries in the html_alignment for all new families
    '''

    query = ("SELECT ha.type, count(*)\n"
             "FROM html_alignment ha, family f\n"
             "WHERE f.rfam_acc=ha.rfam_acc\n"
             "AND f.author like \'%Arga%\'\n"
             "GROUP BY ha.type")

    # TO BE IMPLEMENTED

# -----------------------------------------------------------------------------


def print_report(no_fams):
    '''
        Calls all functions and displays the results on screen
    '''

    cnx = RfamDB.connect()

    check_ss_images(cnx, no_fams)
    check_sunburst(cnx)
    count_rchie_diagrams(cnx, no_fams)
    check_alignment_and_tree(cnx, no_fams)
    check_html_alignment(cnx, no_fams)

    RfamDB.disconnect(cnx)


# -----------------------------------------------------------------------------

def usage():
    '''
        Displays usage information on screen
    '''

    # TO BE IMPLEMENTED

    pass

# -----------------------------------------------------------------------------
if __name__ == '__main__':

    # maybe no_fams not required or provide a number of params to run the test on
    # maybe calculate no_fams from DB
    no_fams = sys.argv[1]
    print_report(no_fams)
