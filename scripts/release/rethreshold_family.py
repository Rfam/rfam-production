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
MEMORY = 2000 
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

	mutually_exclusive_args=parser.add_mutually_exclusive_group()
    	mutually_exclusive_args.add_argument('-f', help='a file containing a list of Rfam family accessions', type=str)
    	mutually_exclusive_args.add_argument('--all', help='runs rfsearch on all families', action="store_true")
    	mutually_exclusive_args.add_argument('--acc', help="a valid rfam family accession RFXXXXX",
                        type=str, default=None)
	
	parser.add_argument('--rfmake', help='run rfmake after rfsearch completion', 
			   action="store_true")
	parser.add_argument('-v', help='runs validation checks', action="store_true")
	
	parser.add_argument('--report', help='generates search reports', action="store_true")

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

def count_hits(scores_file):
	"""
	Function to count SEED and FULL hits in outlist and species files at three
	different thresholds (above ga, below ga, below rev)

	scores_file: This is either the species or the outlist files from the family
	directories

	return: A dictionary with SEED and FULL counts at different thresholds
	"""

	# check point flags
	flag_curr = 0
	flag_rev = 0

	# initialization of counts
	counts = {"seed_above_ga": 0,
		"full_above_ga": 0,
		"full_below_ga": 0,
		"seed_below_ga": 0,
		"seed_below_rev": 0,
		"full_below_rev": 0 }

	# load file for easy parsing
        fp = open(scores_file, 'r')

	# generate stats
	for line in fp:
		# make flag_curr = 1 when we reach that line
		if line.find("CURRENT THRESHOLD") != -1:
			flag_curr = 1
			continue

		# when we reach the reversed sequence line set the flag to 1
		if line.find("BEST REVERSED") != -1:
			flag_rev = 1
			continue

		# we are above the 
		if flag_curr == 0 and flag_rev == 0:
			if line.find("SEED") != -1:
				counts["seed_above_ga"] += 1
		
			elif line.find("FULL") != -1:
				counts["full_above_ga"] += 1
		# we are somewhere in between current threshold and reversed cutoff
		elif flag_curr == 1 and flag_rev == 0:
			if line.find("SEED") != -1:
				counts["seed_below_ga"] += 1
		
			elif line.find("FULL") != -1:
				counts["full_below_ga"] += 1

		elif flag_curr == 1 and flag_rev == 1:
			if line.find("SEED") != -1:
				counts["seed_below_rev"] += 1
		
			elif line.find("FULL") != -1:
				counts["full_below_rev"] += 1
		
	fp.close()
	
	
	return counts

# ----------------------------------------------------------------------------------

def generate_search_stats(family_dir, scores_file = 'species'):
	"""
	Function to generate useful search stats per family

	family_dir: A valid Rfam family checkout directory where pre-computed searches
	were ran
	scores_file: A string specifying the scores file to parse (outlist, species) 

	return: report string
	"""
	
	rfam_acc = os.path.basename(family_dir)
	# check point flags
        flag_curr = 0
        flag_rev = 0
	elements = None
	prev_line = None
        
	# initialization of counts
        counts = {"seed_above_ga": 0,
                "full_above_ga": 0,
                "full_below_ga": 0,
                "seed_below_ga": 0,
                "seed_below_rev": 0,
                "full_below_rev": 0 }

	# get some useful numbers from the database
        no_seed_seqs_db = db.get_number_of_seed_sequences(rfam_acc)
	no_full_hits_db = db.get_number_of_full_hits(rfam_acc)
	unique_ncbi_ids_db = db.get_family_unique_ncbi_ids(rfam_acc)

	scores_fp = open(os.path.join(family_dir, scores_file), 'r')

	# this will basically read the first line which is a header so no harm
	line = scores_fp.readline()
	prev_line = line

	ncbi_ids_from_seed = []
	
	# generate stats
        for line in scores_fp:
		
		# make flag_curr = 1 when we reach that line
                if line.find("CURRENT THRESHOLD") != -1:
                        flag_curr = 1
                       	continue

                # when we reach the reversed sequence line set the flag to 1
                if line.find("BEST REVERSED") != -1:
                        flag_rev = 1
                        continue

		if line[0] != '#':
                        elements = [x for x in line.strip().split(' ') if x!= '']

			# add id to ncbi_ids
			ncbi_ids_from_seed.append(elements[5])
                	
			# we are above the GA
                	if flag_curr == 0 and flag_rev == 0:
                        	if elements[2] == "SEED":
                                	counts["seed_above_ga"] += 1

                        	elif elements[2] == "FULL":
                                	counts["full_above_ga"] += 1
                
			# we are somewhere in between current threshold and reversed cutoff
                	elif flag_curr == 1 and flag_rev == 0:
				if elements[2] == "SEED":
                                	counts["seed_below_ga"] += 1

                        	elif elements[2] == "FULL":
                                	counts["full_below_ga"] += 1

                	elif flag_curr == 1 and flag_rev == 1:
                        	if elements[2] == "SEED":
                                	counts["seed_below_rev"] += 1

                      		elif elements[2] == "FULL":
                                	counts["full_below_rev"] += 1

        scores_fp.close()

	new_ncbi_ids_found = len(list(set(unique_ncbi_ids_db).difference(set(ncbi_ids_from_seed))))

	print ("rfam_acc\tnum_seed_seqs\tnum_full_hits\tnum_ncbi_ids\tnum_seed_above_GA\tnum_full_above_ga\t"),
	print ("tnum_seed_below_ga\tnum_full_below_ga\tSEED_below_rev\tfull_below_rev\tnum_new_ncbi_ids")
	print ('\t'.join([rfam_acc, str(no_seed_seqs_db), str(no_full_hits_db), str(len(unique_ncbi_ids_db)), 
			str(counts["seed_above_ga"]), str(counts["full_above_ga"]), str(counts["seed_below_ga"]), 
			str(counts["full_below_ga"]), str(counts["seed_below_rev"]), str(counts["full_below_rev"]), 
			str(new_ncbi_ids_found)]))

# ----------------------------------------------------------------------------------

def write_family_report_file(family_dir, scores_file = "species"):
	"""
	Function to generate a report about the outcome of a new search 

	family_dir: A valid location of an Rfam family checkout
	scores_file: This is a string which specifies the file to parse (outlist | species)
	It parses species file by default.

	return (int): A number specifying the curation priority for a specific family, where
	3: critical, 2: critical but not erroneous, 1: check seed, 0: no attention needed    
	"""

	priority = 0
	# fetch number of seed sequences from the database
	rfam_acc = os.path.basename(family_dir)
	no_seed_seqs = db.get_number_of_seed_sequences(rfam_acc)	

	scores_file_loc = os.path.join(family_dir, scores_file)
	counts = count_hits(scores_file_loc)

	report_fp = open(os.path.join(family_dir, "search_report.txt"),'w')

	# sum all seed counts to get total number of seed sequences
	counted_seed_seqs = counts["seed_above_ga"] + counts["seed_below_ga"] + counts["seed_below_rev"]
	
	# Critical SEED issues
	if counts["seed_below_rev"] != 0:
		report_fp.write("CRITICAL: %s SEED sequences below reversed cutoff\n" % str(counts["seed_below_rev"]))
		priority = 3

	if counts["seed_below_ga"] > counts["seed_above_ga"]:
		percentage = float(counts["seed_below_ga"] * 100) / float(no_seed_seqs)
		report_fp.write("CRITICAL: More SEED sequences below GA than above. %s\n" % percentage)
		
		if priority < 2:
			priority = 2
	
	if counted_seed_seqs != no_seed_seqs:
		report_fp.write("WARNING: The number of SEED sequences in the database does not match the number in the alignment\n\n")			
		priority = 3

	# TODO - Develop code to check taxonomic distribution
        # TODO - Use information from FULL hits too	

	# some useful information
	report_fp.write("Total number of SEED sequences in DB: %s\n" % no_seed_seqs)
	report_fp.write("Total number of SEED sequences counted: %s\n" % counted_seed_seqs)

	report_fp.write("%s SEED sequences are above GA\n" % counts["seed_above_ga"])
	report_fp.write("%s SEED sequences are below GA\n" % counts["seed_below_ga"])
	report_fp.write("%s SEED sequences are below the reversed cutoff\n" % counts["seed_below_rev"])	
	
	report_fp.close()
	
	return priority

# ----------------------------------------------------------------------------------

if __name__ == '__main__':

    # create a new argument parser object
    parser = parse_arguments()
    args = parser.parse_args()

    family_dir = "/hps/nobackup/production/xfam/rfam/RELEASES/14.2/family_searches/batch1/RF00721"
    generate_search_stats(family_dir, scores_file = 'species')
    sys.exit()

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

    elif args.report:
	if args.acc:
		# check if searches where validated 
		if not os.path.exists(os.path.join(args.dest_dir, "validation.log")):
			sys.exit("WARNING: This search may be invalid. Run validation and try again!")
		family_dir = os.path.join(args.dest_dir, args.acc)
		species_file = os.path.join(family_dir, "species")
	elif args.all:
		families  = [x for x in os.listdir(args.dest_dir) if os.path.isdir(x)]
		for family in families:
			family_dir = os.path.join(args.dest_dir)
			generate_search_stats(family_dir, scores_file = 'species')
		
