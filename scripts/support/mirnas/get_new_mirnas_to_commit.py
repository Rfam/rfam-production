import os
import json

if __name__=='__main__':

	to_commit_fp = "/hps/nobackup/production/xfam/rfam/RELEASES/14.3/miRNA_relabelled/support_code/data/new_mirnas_2020-11-17.json"
	committed_fp = "/hps/nobackup/production/xfam/rfam/RELEASES/14.3/miRNA_relabelled/support_code/data/committed_mirnas.json"
	
	fp = open(to_commit_fp, 'r')
	to_commit = json.load(fp)
	fp.close()	

	fp = open(committed_fp, 'r')
	committed = json.load(fp)
	fp.close()

	novel = {}

	for mirna in to_commit.keys():
		if mirna not in committed:
			if mirna not in novel:
				novel[mirna] = ""
		else:
			print mirna
	print "to_commit: ", len(to_commit.keys())
	print "committed: ", len(committed.keys())
	print "novel: ", len(novel.keys())
