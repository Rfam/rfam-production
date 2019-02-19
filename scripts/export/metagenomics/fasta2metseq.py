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

"""
Script to generate a metseq dump
"""

import os
import sys
import subprocess as p

# ---------------------------------------------------------------------------------------------------------------

def fasta_to_metseq_dump(fasta_input, filename=None, dest_dir=None):
	"""
	Convert UMGS fasta file header accessions to metseq txt dump
	
	fasta_input: A UMGS fasta file or a directory with multiple fasta files
	filename: The output filename. If None, uses the UMGS filename by default
	dest_dir: The destination directory. If None, uses the input directory by default

	returns: Void
	"""

	if dest_dir is None:
		dest_dir =os.path.split(fasta_input)[0]

	if os.path.isdir(fasta_input):
		
		# make sure we only list the fasta files
		fasta_files = [x for x in os.listdir(fasta_input) if x.endswith(".fa") or x.endswith(".fasta")]

		for fasta_file in fasta_files:
			filename = fasta_file.partition('.')[0]

			fasta_file_loc = os.path.join(fasta_input, fasta_file)
			temp_file_loc = os.path.join(fasta_input, "."+filename+".accs")
			# create a hidden accession file
			p.call("grep \">\" %s > %s" % (fasta_file_loc, temp_file_loc), shell=True)
			
			fp_out = open(os.path.join(dest_dir, filename+'.metseq'), 'w')
			fp_in = open(temp_file_loc, 'r')

			for line in fp_in:
				# strip whitespaces and trim off '>' char
				metseq_acc = line.strip()[1:]
				fp_out.write(filename+'\t'+metseq_acc+'\n')

			fp_out.close()
			fp_in.close()

			# remove temp file
			os.remove(temp_file_loc)

	elif os.path.isfile(fasta_input):
		# get filename
		input_components = os.path.split(fasta_input)
		
		if filename is None:	
                	filename = input_components[1].partition('.')[0]
		
		temp_file_loc = os.path.join(input_components[0], "."+filename+".accs")
		
		p.call("grep \">\" %s > %s" % (fasta_input, temp_file_loc), shell=True)

		fp_out = open(os.path.join(dest_dir, filename+'.metseq'), 'w')
		fp_in = open(temp_file_loc, 'r')

		for line in fp_in:
			metseq_acc = line.strip()[1:]
			fp_out.write(filename+'\t'+metseq_acc+'\n')

		fp_out.close()
		fp_in.close()
	
		os.remove(temp_file_loc)
		

# ---------------------------------------------------------------------------------------------------------------

if __name__ == '__main__':
	
	fasta_input = sys.argv[1]
	dest_dir = sys.argv[2]

	fasta_to_metseq_dump(fasta_input, None, dest_dir)
