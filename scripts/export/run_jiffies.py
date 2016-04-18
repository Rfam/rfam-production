#!/usr/bin/python

'''
Created on 29 Feb 2016

@author: ikalvari

Description: Support code designed to ease Rfam jiffies' execution.

Notes: 1. Call this script as rfamprod
'''

# ---------------------------------IMPORTS-------------------------------------

import os
import sys
import subprocess

# -----------------------------------------------------------------------------


def call_jiffy(jiffy, fam_file, outdir=None):
    '''
        This function was designed to call the perl script defined by the jiffy 
        parameter. Used in Rfam 12.1 with jiffies writeAnnotatedCM.pl and
        writeAnnotatedSeed.pl to generate CM and SEED files for the FTP server.

        jiffy: The path to the perl script to call
        fam_file: A list of all rfam_accessions
        outdir: Destination directory where the files will be generated

    '''

    # TO DO
    # convert this to load family accessions from the database

    # change directory to outdir. Due to perl script limitations, to have the
    # output generated in the specified directory, the jiffy needs to be called
    # inside that
    if outdir is not None:
        os.chdir(outdir)

    fp = open(os.path.abspath(fam_file), 'r')

    for rfam_acc in fp:
        cmd = "%s %s" % (jiffy, rfam_acc)
        subprocess.call(cmd, shell=True)

    fp.close()

# -----------------------------------------------------------------------------


def usage():

        # TO DO
    pass
# -----------------------------------------------------------------------------
if __name__ == '__main__':

    # need to check the number of parameters provided
    jiffy = sys.argv[1]
    fam_file = sys.argv[2]
    outdir = sys.argv[3]

    call_jiffy(jiffy, fam_file, outdir)
