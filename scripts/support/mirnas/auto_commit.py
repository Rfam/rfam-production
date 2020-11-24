import os
import sys
import subprocess
import json

from subprocess import Popen, PIPE

search_dirs = ["/hps/nobackup/production/xfam/rfam/RELEASES/14.3/miRNA_relabelled/batch1_chunk1_searches",
              "/hps/nobackup/production/xfam/rfam/RELEASES/14.3/miRNA_relabelled/batch1_chunk2_searches",
              "/hps/nobackup/production/xfam/rfam/RELEASES/14.3/miRNA_relabelled/batch2/searches"]

# ---------------------------------------------------------------------------------------------


def check_desc_ga(DESC, cut_ga):
	"""
	"""

	process = Popen(['grep', "GA", DESC], stdin=PIPE, stdout=PIPE, stderr=PIPE)
        output, err = process.communicate()
	
        if output.find("%.2f"%float(cut_ga)) == -1:
                return False

        return True 
	

# ---------------------------------------------------------------------------------------------


def check_family_passes_qc(family_dir):
	
	dir_elements = os.path.split(family_dir)
	search_dir = dir_elements[0]

	os.chdir(search_dir)

	process = Popen(["rqc-all.pl", dir_elements[1]], stdin=PIPE, stdout=PIPE, stderr=PIPE)
        output = process.communicate()[1]

    	if output.find("Family passed with no serious errors") == -1:
        	return False

    	return True 

# ---------------------------------------------------------------------------------------------


def commit_family(family_dir, mirna_name):

	dir_elements = os.path.split(family_dir)
	os.chdir(dir_elements[0])
	family_dir = dir_elements[1]

	process = Popen(['rfnew.pl', '-m', "\"Adding new miRNA family %s \""% (mirna_name), family_dir], stdin=PIPE, stdout=PIPE, stderr=PIPE)
    	output = process.communicate()[1]

    	if output.find("This family has been asssigned the accession") == -1:
        	return False

    	return True

# ---------------------------------------------------------------------------------------------


if __name__=='__main__':

	fp = open(sys.argv[1], 'r')	
	miRNA_accessions = json.load(fp)
	fp.close()

	fp = open("/hps/nobackup/production/xfam/rfam/RELEASES/14.3/input/existing_mirna_families.json", 'r')
	existing_fams = json.load(fp)
	fp.close()

	committed = {}
	
	skip = ["MIPF0001496__mir-6012", "MIPF0001508__mir-4504", "MIPF0001511__mir-4803"]
	
	for miRNA in skip:
		del(miRNA_accessions[miRNA])
	
	fp = open(os.path.join(os.getcwd(), 'successful_mirna_commits.log'), 'w')

	for accession in miRNA_accessions.keys():
		if accession not in committed:
			dir_label = ''
			if accession.find("_relabelled")==-1:
				dir_label = accession+"_relabelled"
			
			for search_dir in search_dirs:
				family_dir_loc = os.path.join(search_dir, dir_label)
				if os.path.exists(family_dir_loc):
					desc_file = os.path.join(family_dir_loc, "DESC")
					if check_desc_ga(desc_file, miRNA_accessions[accession]) is True:
						if check_family_passes_qc(family_dir_loc) is True:
							mirna_name = ""
							if accession[0:2]=='MI':
								mirna_name = accession.split("_")[2]
							else:
								mirna_name = accession.split("_")[0]
							check = commit_family(family_dir_loc, mirna_name)
							if check is True:
								committed[accession] = ""
							print "Family %s committed" % (accession)
							fp.write(accession+'\n')
				else:
					continue
	fp.close()
	fp = open("/hps/nobackup/production/xfam/rfam/RELEASES/14.3/miRNA_relabelled/support_code/committed_mirnas.json", 'w')
	json.dump(committed, fp)
	fp.close()
