"""
Copyright [2009-2019] EMBL-European Bioinformatics Institute
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

# ----------------------------------------------------------------------------------

import os
import sys
import subprocess
import argparse

# ------------------------------------- GLOBALS ------------------------------------

LSF_GROUP = "/family_srch"
MEMORY = 8000
CPU = 8

# ----------------------------------------------------------------------------------


def checkout_family(rfam_acc):
    """
    Checks out a family from Rfam based on a valid Rfam accession.

    rfam_acc: A valid Rfam accession

    return: None
    """

    cmd = "rfco.pl %s" % rfam_acc

    subprocess.call(cmd, shell=True)

    # add some checks here

# ----------------------------------------------------------------------------------


def submit_new_rfsearch_job(family_dir):
    """
    Submits a new lsf job that runs rfsearch to update SCORES for a new release.
    If no threshold is set with rfsearch.pl, it uses existing thresholds by default.

    family_dir: The physical location of the family directory

    return: None
    """
    # use the pre-process command to change directory to family_dir

    rfam_acc = os.path.basename(family_dir)

    lsf_err_file = os.path.join(family_dir, "auto_rfsearch.err")
    lsf_out_file = os.path.join(family_dir, "auto_rfsearch.out")

    cmd = ("bsub -M %s -R \"rusage[mem=%s]\" -o %s -e %s -n %s -g %s -q production-rh7 "
          "-J %s cd %s && rfsearch.pl -cnompi")

    subprocess.call(cmd % (MEMORY, MEMORY, lsf_out_file, lsf_err_file,
                         CPU, LSF_GROUP, rfam_acc, family_dir), shell=True)

# ----------------------------------------------------------------------------------


def load_rfam_accessions_from_file(accession_list):
    """
    This function parses a .txt file containing Rfam accessions and returns those
    accession_list: This is a .txt file containing a list of Rfam accessions

    return: list of Rfam family accessions
    """
    
    fp = open(accession_list, 'r')
    
    accessions = [x.strip() for x in fp]

    fp.close()

    return accessions

# ----------------------------------------------------------------------------------


def checkout_and_search_family(rfam_acc, dest_dir):
	"""
	This function combines family checkout (rfco.pl) and re-scoring of hits 
	using rfsearch.pl. 

	rfam_acc: A valid Rfam family accession (RFXXXXX)
	dest_dir: A valid destination directory, where to checkout the family

	return: void
	"""

	# get family directory
	family_dir = os.path.join(dest_dir, rfam_acc)
	# checkout family if not done already
        if not os.path.exists(family_dir):
        	checkout_family(rfam_acc)
	
	submit_new_rfsearch_job(family_dir)

# ----------------------------------------------------------------------------------


def parse_arguments():
	"""
	Uses python's argparse to parse the command line arguments

	return: Argparse parser object
	"""

	# create a new argument parser object
    	parser = argparse.ArgumentParser(description='Update scores for new release')

    	# group required arguments together
    	req_args = parser.add_argument_group("required arguments")
    	req_args.add_argument('--dest_dir', help='destination directory where to checkout families',
                        type=str, required=True)

    	parser.add_argument('-f', help='a file containing a list of Rfam family accessions', type=str)
    	parser.add_argument('--all', help='runs rfsearch on all families', type=str)
    	parser.add_argument('--acc', help="a valid rfam family accession RFXXXXX",
                        type=str, default=None)

	return parser

# ----------------------------------------------------------------------------------


if __name__ == '__main__':

    # create a new argument parser object
    parser = parse_arguments()
    args = parser.parse_args()

    #if args.h:
    #    parser.print_help()

    if args.acc:
        if args.acc[0:2] == 'RF' and len(args.acc) == 7:
            os.chdir(args.dest_dir)
   	    checkout_and_search_family(args.acc, args.dest_dir)         
    
    elif args.f:
	if not os.path.isfile(args.f):
		print "The file location you provided does not exist!\n"
	# move to destination directory
	os.chdir(args.dest_dir)
	accessions = load_rfam_accessions_from_file(args.f)

	for rfam_acc in accessions:
		checkout_and_search_family(rfam_acc, args.dest_dir)		
	
