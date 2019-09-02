import os
import sys
import json
import subprocess
from config import rfam_config as gf

# -----------------------------------------------------------------------------
# TO DO, set upid_taxid_file to None and try to fetch ncbi ids from Uniprot.
# There's a chance that upids become unavailable though
def extract_metadata_from_fasta(fasta_file, taxid, source, filename=None, to_file=True):
    """
    Parses a fasta file and generates rfamseq like matadata using esl-seqstat

    fasta_file: A valid fasta file
    taxid: A file with upid and taxid mappings.
    database:

    returns: void
    """

    mol_type = "genomic DNA"
    previous_acc = ''

    # create an output file pointed
    output_fp = None

    if to_file is True:
        if filename is None:
            filename = os.path.basename(fasta_file).partition('.')[0]
        destination = os.path.split(fasta_file)[0]
        output_fp = open(os.path.join(destination, filename+".rfamseq"), 'w')

    args = [gf.ESL_SEQSTAT, "-a", "--dna", fasta_file]
    # open a process pipe and run esl-seqstat. Stores the result in process
    process = subprocess.Popen(args, stdout=subprocess.PIPE)
    
    # fetch esl-seqstat results
    result = process.communicate()

    # read in and clean-up sequence stats
    seq_stats = [x[2:] for x in result[0].strip().split('\n') if x[0] == '=']

    # process every sequence header to generate a rfamseq entry
    for sequence in seq_stats:
        seq_elements = [x.strip() for x in sequence.split(' ') if x != '']

        rfamseq_acc = seq_elements[0].split('|')[-1]
        accession = rfamseq_acc.partition('.')[0]
        version = rfamseq_acc.partition('.')[2]
        length = seq_elements[1]
        description = " ".join(seq_elements[2:])

        rfamseq_entry = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (rfamseq_acc, accession,
                                                      version, taxid,
                                                      mol_type, length,
                                                      description, previous_acc,
                                                      source)
        if to_file is True:
            output_fp.write(rfamseq_entry+'\n')
        else:
            print rfamseq_entry

    if to_file is True:
        output_fp.close()

# -----------------------------------------------------------------------------


def main(fasta_input, upid_list, upid_gca_tax_file):
    """
    This is the main function of the fasta2rfamseq script, which converts a
    fasta file to rfamseq table dumps (e.g. filename.rfamseq) for easy import
    to rfam_live upon release

    fasta_input: This can be a single fasta file or a genome project directory
	         as orgnised by the genome download pipeline
    upid_list:   This is a plain txt file listing all the upids in the genome
                 project directory for which to generate the rfamseq files. 
                 It defaults to None, for single fasta files
    upid_gca_tax_file: This a json file mapping each upid to its corresponding
                       taxid. In the case of Zasha or RNAcentral input, this 
		       will be a mapping of the sequence accession in the 
  		       fasta file to its associated taxid
    
    return: void
    """

    fp = open(upid_gca_tax_file, 'r')
    upid_gca_tax_dict = json.load(fp)
    fp.close()

    if os.path.isfile(upid_list):
    	project_dir = fasta_input
	fp = open(upid_list, 'r')
    	upids = [x.strip() for x in fp]
    	fp.close()

    	for upid in upids:
        	print upid
        	subdir = os.path.join(project_dir, upid[-3:])
        	updir = os.path.join(subdir, upid)
        	upfasta = os.path.join(updir, upid+'.fa')

        	# do some sanity checks
        	if os.path.exists(upfasta):
            	# defining the source of the sequences according to assembly accession
            		source = "UNIPROT; ENA"
            		if upid in upid_gca_tax_dict:
	    			if upid_gca_tax_dict[upid]["GCA"] != -1:
                			if upid_gca_tax_dict[upid]["GCA"][0:3] == "GCF":
                    				source = "NIH; NCBI"
                			else:
                    				source = "UNIPROT; ENA"

	    	else:
			print "%s not in the current Uniprot release"%upid
			continue

	    	extract_metadata_from_fasta(upfasta, upid_gca_tax_dict[upid]["TAX"],
                                        source, filename=None, to_file=True)

    else:
    	project_dir = fasta_input
    	upid = upid_list
    	subdir = os.path.join(project_dir, upid[-3:])
            updir = os.path.join(subdir, upid)
            upfasta = os.path.join(updir, upid+'.fa')

    	source = "UNIPROT; ENA"
	if upid in upid_gca_tax_dict:
        	if upid_gca_tax_dict[upid]["GCA"] != -1:
			if upid_gca_tax_dict[upid]["GCA"][0:3] == "GCF":
				source = "NIH; NCBI"
			else:
				source = "UNIPROT; ENA"
	else:
		print "%s not in the current Uniprot release" % upid
		sys.exit()

	extract_metadata_from_fasta(upfasta, upid_gca_tax_dict[upid]["TAX"],
                                        source, filename=None, to_file=True)


# -----------------------------------------------------------------------------


if __name__ == '__main__':

    project_dir = sys.argv[1]
    upid_list = sys.argv[2]
    upid_gca_tax_file = sys.argv[3]

    main(project_dir, upid_list, upid_gca_tax_file)
