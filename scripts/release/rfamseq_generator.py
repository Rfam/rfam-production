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

import os
import sys
import re
import argparse

# ------------------------------------------------------------------------


def merge_all_genome_files(project_dir, dest_dir, filename='rfamseq'):
    """
    Simple script to merge all genomes to a single rfamseq file

    project_dir: The path to a genome download project directory
    dest_dir: The directory where to create the new rfamseq file
    filename: A filename for the rfamseq file. Defaults to rfamseq

    return: Void
    """

    err_cases_fp = os.path.join(dest_dir, filename+'_err_cases.txt')

    rfamseq_fp = open(os.path.join(dest_dir, filename + ".fa"), 'w')

    subdirs = [x for x in os.listdir(project_dir)
               if os.path.isdir(os.path.join(project_dir, x))]

    for subdir in subdirs:
        subdir_loc = os.path.join(project_dir, subdir)
        updirs = os.listdir(subdir_loc)

        for upid in updirs:
            updir_loc = os.path.join(subdir_loc, upid)
            up_fasta = os.path.join(updir_loc, upid + ".fa")

            if os.path.exists(up_fasta):
                fasta_fp = open(up_fasta, 'r')
                for seq_line in fasta_fp:
                    #if seq_line[0] == '>':
                    rfamseq_fp.write(seq_line)
                    #else:
                    #    if seq_validator(seq_line):
                    #        rfamseq_fp.write(seq_line)
                    #    else:
                    #        print upid + '\t' + seq_line

                fasta_fp.close()

    rfamseq_fp.close()

# ------------------------------------------------------------------------

def merge_project_files(project_dir, dest_dir, file_type, filename='rfamseq'):
    """
    Simple script to merge all genomes to a single rfamseq file

    project_dir: The path to a genome download project directory
    dest_dir: The directory where to create the new rfamseq file
    filename: A filename for the rfamseq file. Defaults to rfamseq

    return: Void
    """

    err_cases_fp = os.path.join(dest_dir, filename+'_err_cases.txt')

    rfamseq_fp = open(os.path.join(dest_dir, filename + "." + file_type), 'w')

    subdirs = [x for x in os.listdir(project_dir)
               if os.path.isdir(os.path.join(project_dir, x))]

    for subdir in subdirs:
        subdir_loc = os.path.join(project_dir, subdir)
        updirs = os.listdir(subdir_loc)

        for upid in updirs:
            updir_loc = os.path.join(subdir_loc, upid)
            up_fasta = os.path.join(updir_loc, upid + "." + file_type)

            if os.path.exists(up_fasta):
                fasta_fp = open(up_fasta, 'r')
                for seq_line in fasta_fp:
                    #if seq_line[0] == '>':
                    rfamseq_fp.write(seq_line)
                    #else:
                    #    if seq_validator(seq_line):
                    #        rfamseq_fp.write(seq_line)
                    #    else:
                    #        print upid + '\t' + seq_line

                fasta_fp.close()

    rfamseq_fp.close()

# ------------------------------------------------------------------------


def seq_validator(sequence):
    """
    Checks if the sequence provided is valid fasta sequence. Returns True
    if the sequence is valid, otherwise returns False.

    sequence: A string for validation
    """

    # checks for ascii characters that should not appear in a fasta sequence
    seq_val = re.compile("[^ATKMBVCNSWD-GUYRHatkbbvcnswdguyrh]")

    # if any illegal characters found return False
    if seq_val.search(sequence):
        return False

    return True

# ------------------------------------------------------------------------

def merge_files_from_accession_list(project_dir, acc_list_file, dest_dir, file_type, filename='rfamseq'):
    """
    Simple script to merge all genomes to a single rfamseq file

    project_dir: The path to a genome download project directory
    dest_dir: The directory where to create the new rfamseq file
    filename: A filename for the rfamseq file. Defaults to rfamseq

    return: Void
    """
    
    if file_type.lower() == 'fasta':
	file_type = 'fa'
    elif file_type.lower() == 'tblout':
	file_type = 'tbl'

    err_cases_fp = os.path.join(dest_dir, filename+'_err_cases.txt')

    rfamseq_fp = open(os.path.join(dest_dir, filename + "." + file_type), 'w')

    #subdirs = [x for x in os.listdir(project_dir)
    #           if os.path.isdir(os.path.join(project_dir, x))]

    fp_in = open(acc_list_file, 'r')
    upids = [x.strip() for x in fp_in]
    fp_in.close()

    #for subdir in subdirs:
    #    subdir_loc = os.path.join(project_dir, subdir)
    #    updirs = os.listdir(subdir_loc)

    for upid in upids:
	suffix = upid[-3:]
	subdir_loc = os.path.join(project_dir, suffix)
	
	if os.path.exists(subdir_loc):
		updir_loc = os.path.join(subdir_loc, upid)

		if os.path.exists(updir_loc):
			up_fasta = os.path.join(updir_loc, upid + "." + file_type)

            		if os.path.exists(up_fasta):
                		fasta_fp = open(up_fasta, 'r')
                		for seq_line in fasta_fp:
                    		#if seq_line[0] == '>':
                    			rfamseq_fp.write(seq_line)
                    		#else:
                    		#    if seq_validator(seq_line):
                    		#        rfamseq_fp.write(seq_line)
                    		#    else:
                    		#        print upid + '\t' + seq_line

                		fasta_fp.close()
    rfamseq_fp.close()

# ------------------------------------------------------------------------

def parse_arguments():
	"""
	Uses python's argparse to parse the command line arguments
	
	return: Argparse parser object
	"""

	# create a new argument parser object
    	parser = argparse.ArgumentParser(description='Merges genomes into a unified fasta file (Rfamseq)')

    	# group required arguments together
    	req_args = parser.add_argument_group("required arguments")
    	req_args.add_argument('--project_dir', help='a project directory where the genome directories reside',
                        type=str, required=True)
	
    	parser.add_argument('-f', help='a file containing a list of valid UPIDs', type=str)
    	parser.add_argument('--all', help='merges all genomes to build the new Rfamseq', action="store_true")
	parser.add_argument('--out', help='a name for the output and or full destination path', type=str)
	parser.add_argument('--type', help='the type of the file to merge (fasta, tblout, rfamseq, genseq)', 
				choices=['fasta', 'tblout', 'rfamseq', 'genseq'], type=str, default='fasta')
	
	return parser

# ------------------------------------------------------------------------
	
if __name__ == '__main__':

    # create a new argument parser object
    parser = parse_arguments()
    args = parser.parse_args()

    dest_dir = os.getcwd()
    filename = 'rfamseq.fasta'

    if os.path.isfile(args.out) or args.out.find("/")!= -1 :
	absolute_path = os.path.abspath(args.out)
	info = os.path.split(args.out)
	dest_dir = info[0]
	filename = info[1]
    else:
	filename = args.out

    if args.all is True:
	project_dir = args.project_dir
	file_type = args.type

    	#merge_all_genome_files(project_dir, dest_dir, filename=filename)
    	merge_project_files(project_dir, dest_dir, file_type, filename)
    
    elif args.f and args.all is False:
	project_dir = args.project_dir
	acc_list_file = os.path.abspath(args.f)
        file_type = args.type

	merge_files_from_accession_list(project_dir, acc_list_file, dest_dir, file_type, filename)
