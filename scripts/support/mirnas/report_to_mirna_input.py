import os
import sys
import json
from datetime import date

import argparse

# -------------------------------------------------------------------------------


def extract_new_mirnas_from_report(report_tsv, type='new'):
	"""
	"""

	new_mirnas = {}

	fp = open(report_tsv, 'r')
        count = 0	
	for line in fp:
		line = line.strip().split('\t')
	        # check if candidate mirna is a new family
		if line[6].lower() == "new family":
			# skip families requiring review
			if line[1] != '?' and line[1]!= '' and line[2]!="1_SEED":		
				if line[0] not in new_mirnas:
					print line
					new_mirnas[line[0]] = line[1]
		elif line[6].lower() == 'done':
			count+=1
	fp.close()
	print count
	return new_mirnas


# -------------------------------------------------------------------------------


def parse_arguments():
	"""
	"""

	parser = argparse.ArgumentParser()

	parser.add_argument("--report", help="miRNA report in .tsv format", action="store")
	parser.add_argument("--dest-dir", help="Desctination directory", action="store", default=os.getcwd())
	
	return parser


# -------------------------------------------------------------------------------


if __name__=='__main__':

	parser = parse_arguments()
	args = parser.parse_args()	

	new_mirnas = extract_new_mirnas_from_report(args.report, type='new')

	fp_out = open(os.path.join(args.dest_dir, "new_mirnas_"+str(date.today())+".json"), 'w')

	json.dump(new_mirnas, fp_out)
	
	fp_out.close()

