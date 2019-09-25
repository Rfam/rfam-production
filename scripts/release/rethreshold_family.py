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
import math
import subprocess
import argparse
from subprocess import Popen, PIPE

from utils import db_utils as db
# ------------------------------------- GLOBALS ------------------------------------
# this group only allows 10 rfsearch jobs to run concurrently
# this means 10*100 = 1000 jobs running concurrently which is the lsf limit
 
LSF_GROUP = "/family_srch"
MEMORY = 8000
CPU = 8
MAX_JOB_COUNT = 1000

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


def submit_new_rfsearch_job(family_dir, rfmake=False):
    """
    Submits a new lsf job that runs rfsearch to update SCORES for a new release.
    If no threshold is set with rfsearch.pl, it uses existing thresholds by default.

    family_dir: The physical location of the family directory
    rfmake: If True, run rfmake after rfsearch completes. Default False
 
    return: None
    """
    # use the pre-process command to change directory to family_dir

    rfam_acc = os.path.basename(family_dir)

    lsf_err_file = os.path.join(family_dir, "auto_rfsearch.err")
    lsf_out_file = os.path.join(family_dir, "auto_rfsearch.out")

    cmd = ("bsub -M %s -R \"rusage[mem=%s]\" -o %s -e %s -n %s -g %s -q production-rh7 "
          "-J %s \"cd %s && rfsearch.pl -cnompi -q production-rh7\"")

    # If rfmake is set to True, runs rfmake following rfsearch, otherwise run rfsearch
    # only by default
    if rfmake is True:
	cmd = ("bsub -M %s -R \"rusage[mem=%s]\" -o %s -e %s -n %s -g %s -q production-rh7 "
          "-J %s \"cd %s && rfsearch.pl -cnompi -q production-rh7 && rfmake.pl\"")

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


def checkout_and_search_family(rfam_acc, dest_dir, rfmake=False):
	"""
	This function combines family checkout (rfco.pl) and re-scoring of hits 
	using rfsearch.pl. If the family directory already exists, then the 
	checkout step will be ignored 

	rfam_acc: A valid Rfam family accession (RFXXXXX)
	dest_dir: A valid destination directory, where to checkout the family
	rfmake: If True, run rfmake after rfsearch completes. Default False 

	return: void
	"""

	# get family directory
	family_dir = os.path.join(dest_dir, rfam_acc)
	# checkout family if not done already
        if not os.path.exists(family_dir):
		os.chdir(dest_dir)
        	checkout_family(rfam_acc)
	
	submit_new_rfsearch_job(family_dir, rfmake)

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
    	parser.add_argument('--all', help='runs rfsearch on all families', action="store_true")
    	parser.add_argument('--acc', help="a valid rfam family accession RFXXXXX",
                        type=str, default=None)
	parser.add_argument('--rfmake', help='run rfmake after rfsearch completion', 
			   action="store_true")
	parser.add_argument('-v', help='runs validation checks', action="store_true")

	return parser

# ----------------------------------------------------------------------------------

def is_valid_family(dest_dir, rfam_acc):
	"""
	Checks if the job ran successfully by checking if .err file is empty and 
	that Success keyword exists in .out file. As an additional sanity check, we
	look for the rfsearch.log file as an indication that rfsearch actually ran.
	
	return: True if the family is valid, False otherwise
	"""

	family_dir = os.path.join(dest_dir, rfam_acc)
	
	# If log file does not exist rfsearch did not run for some reason
	if not os.path.exists(os.path.join(family_dir, "rfsearch.log")):
		return False
	
	# check if lsf .err file is empty 
	if not os.path.getsize(os.path.join(family_dir, "auto_rfsearch.err")) == 0:
		return check_rfsearch_log_success(family_dir)
		#return False

	# check if success in .out file
	lsf_out_file = os.path.join(family_dir, "auto_rfsearch.out")

	process = Popen(['grep', 'Success', lsf_out_file], stdin=PIPE, stdout=PIPE, stderr=PIPE)
	output, err = process.communicate()

	if output.find("Successfully completed.") == -1:
		return False
	
	return True

# ----------------------------------------------------------------------------------

def check_rfsearch_log_success(family_dir):
	"""
	Checks if the rfsearch.log file contains the success string # [ok] in 
	order to mark the family as successfully completed.
	"""

	rfsearch_log_file = os.path.join(family_dir, "rfsearch.log")
	process = Popen(['tail', '-1', rfsearch_log_file], stdin=PIPE, stdout=PIPE, stderr=PIPE)
        output, err = process.communicate()

	if output.find("# [ok]") == -1:
		return False
	
	return True	
# ----------------------------------------------------------------------------------

if __name__ == '__main__':

    # create a new argument parser object
    parser = parse_arguments()
    args = parser.parse_args()

    #if args.h:
    #    parser.print_help()

    if args.acc and not args.v:
	# check accession provided is valid 
        if args.acc[0:2] == 'RF' and len(args.acc) == 7:
            os.chdir(args.dest_dir)
	    
   	    checkout_and_search_family(args.acc, args.dest_dir, rfmake=args.rfmake)
    
    elif args.f and not args.v:
	if not os.path.isfile(args.f):
		print "The file location you provided does not exist!\n"
		sys.exit()
	# move to destination directory
	os.chdir(args.dest_dir)
	accessions = load_rfam_accessions_from_file(args.f)
	
	"""
	# get number of job batches we need to submit
        # casting to int chops off decimals and ceil rounds up to nearest int
	if len(accessions) > MAX_JOB_COUNT:
		no_batches = int(math.ceil(len(accessions)/MAX_JOB_COUNT))
	
	i = 0
	while i < no_batches:
		lidx = i * MAX_JOB_COUNT     # left index
		ridx = (i+1) * MAX_JOB_COUNT # right index
		
		# get exactly MAX_JOB_COUNT items
		if i < no_batches - 1:
			new_batch = accessions[lidx:ridx]
		# get remaining accessions for last batch
		else:
			new_batch = accessions[lidx:]		

		# call function to submit batch
		# while monitoring is True:
		# cluster monitoring function to be called here	
		i+1 # this is done when the monitoring loop becomes false which is a signal to submit another batch
	"""
	for rfam_acc in accessions:
		checkout_and_search_family(rfam_acc, args.dest_dir, rfmake=args.rfmake)		
    
    # run rfsearch on all families in the database
    elif args.all and not args.v:
	# fetch Rfam family accessions from the database
	# call checkout_and_search_family for every family in the list
	# fetches all rfam accessions from the database in DESC order based on the number of sequences in SEEDs
	rfam_acc_list = db.fetch_rfam_accs_sorted(order='DESC')
	for rfam_acc in rfam_acc_list:
		checkout_and_search_family(rfam_acc, args.dest_dir, rfmake=args.rfmake)

    elif args.v:
	if args.acc:
		if not is_valid_family(args.dest_dir, args.acc):
			print "The family %s does not validate!" % args.acc
	elif args.f:
		validation_file = os.path.join(args.dest_dir, "validation.log")
		fp = open(validation_file, 'w')
		accessions = load_rfam_accessions_from_file(args.f)
		for rfam_acc in accessions:
			if not is_valid_family(args.dest_dir, rfam_acc):
				fp.write(rfam_acc + '\n')
		
		fp.close()

		if os.path.getsize(validation_file) == 0:
                        print "Validation process completed! All searches completed successfully!"
                else:
                        print "Validation process completed! Check validation.log for erroneous searches!"

	elif args.all:
		validation_file = os.path.join(args.dest_dir, "validation.log")
		fp = open(validation_file, 'w')
		accessions = [x for x in os.listdir(args.dest_dir) if os.path.isdir(os.path.join(args.dest_dir, x))]
		for rfam_acc in accessions:
			if not is_valid_family(args.dest_dir, rfam_acc):
				fp.write(rfam_acc + '\n')
		
		fp.close()

		if os.path.getsize(validation_file) == 0:
			print "Validation process completed! All searches completed successfully!"
		else:
			print "Validation process completed! Check validation.log for erroneous searches!"
