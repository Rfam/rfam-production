import os
import sys
import subprocess
import argparse

FAMILY_EXCEPTIONS = {'RF02924': '', 'RF03064': '', 'RF02913': '',
                   'RF02543': '', 'RF00017': '', 'RF02540': ''}

# -----------------------------------------------------------------------

def is_outlist_weird(input):

	"""
	"""
	
	if os.path.isdir(input):
 
		families = [x for x in os.listdir(input) if os.path.isdir(os.path.join(input,x))]
	
		for family in families:
			
			if family not in FAMILY_EXCEPTIONS:
				family_dir = os.path.join(input, family)
				outlist_loc = os.path.join(family_dir, "outlist")
				
				if os.path.exists(outlist_loc):
					cmd = "grep REV %s | wc -l" % outlist_loc
					num_rev_lines = int (subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT))

					if num_rev_lines > 1: 
						print family
				else:
					print "ERROR: outlist file for family %s does not exist!" % family

	elif os.path.isfile(input):
		cmd = "grep REV %s | wc -l" % input
		num_rev_lines = int (subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT))

		if num_rev_lines > 1:
			print "ERROR: REVERSED threshold lines > 1. Check outlist!\n"

# -----------------------------------------------------------------------

def parse_arguments():

	"""
	Parse comamnd line arguments
	"""

	parser = argparse.ArgumentParser(description='Checks outlist format')

	arguments = parser.add_argument('--input', help='A directory with multiple family directories or a single outlist file'
					, type=str, required=True)

	return parser

# ----------------------------------------------------------------------- 

if __name__=='__main__':

	args = parse_arguments().parse_args()
	input = args.input
 
	is_outlist_weird(input) 
