import os
import sys
import subprocess

FAMILY_EXCEPTIONS = {'RF02924': '', 'RF03064': '', 'RF02913': '',
                   'RF02543': '', 'RF00017': '', 'RF02540': ''}

# -----------------------------------------------------------------------

def is_outlist_weird(dest_dir):

	"""
	"""

	families = [x for x in os.listdir(dest_dir) if os.path.isdir(os.path.join(dest_dir,x))]
	
	for family in families:
		
		if family not in FAMILY_EXCEPTIONS:
			family_dir = os.path.join(dest_dir, family)
			outlist_loc = os.path.join(family_dir, "outlist")

			cmd = "grep REV %s | wc -l" % outlist_loc
			num_rev_lines = int (subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT))

			if num_rev_lines > 1: 
				print family

# -----------------------------------------------------------------------	

if __name__=='__main__':

	dest_dir = sys.argv[1]
	is_outlist_weird(dest_dir) 
