import os
import sys
import json
import subprocess

searchdirs = ["/hps/nobackup/production/xfam/rfam/RELEASES/14.3/miRNA_relabelled/batch1_chunk1_searches", 
              "/hps/nobackup/production/xfam/rfam/RELEASES/14.3/miRNA_relabelled/batch1_chunk2_searches",
              "/hps/nobackup/production/xfam/rfam/RELEASES/14.3/miRNA_relabelled/batch2/searches"]

MEMORY = 15000
CPU = 4
LSF_GROUP = "/family_srch"


def autorfmake(thresholds_file):
    
    cmd = ("bsub -M %s -o %s -e %s -n %s -g %s -q production-rh74 "
               "-J %s \"cd %s && rfmake.pl -t %s -forcethr -a\"")


    fp = open(thresholds_file, 'r')
    thresholds = json.load(fp)

    for family in thresholds.keys():
        #print family
        for searchdir in searchdirs:
            family_dir = ""

            if family.find("relabelled") == -1:
                family_dir = os.path.join(searchdir, family + "_relabelled")
            else:
                family_dir = os.path.join(searchdir, family)

            if os.path.exists(family_dir):
                lsf_err_file = os.path.join(family_dir, "auto_rfmake.err")
                lsf_out_file = os.path.join(family_dir, "auto_rfmake.out")
		job_name = family

		print family_dir
		#print cmd % (MEMORY, lsf_out_file, lsf_err_file, CPU, LSF_GROUP, job_name, family_dir, str(thresholds[family]))
  		subprocess.call(cmd % (MEMORY, lsf_out_file, lsf_err_file, CPU, LSF_GROUP, job_name, family_dir, str(thresholds[family])), shell=True)
            else:
                continue
	#sys.exit()

# ------------------------------------------------------------------------


def list_mirnas(mirna_dir):
	
	mirnas = []	
	for mirna in os.listdir(mirna_dir):
		if os.path.isdir(os.path.join(mirna_dir, mirna)):
			mirna_label = mirna.replace("_relabelled", "")
			if mirna_label[0]=='M':
				mirnas.append(mirna_label.partition('__')[2])
			else:
				mirnas.append(mirna_label)

	return mirnas

# ------------------------------------------------------------------------


if __name__=='__main__':
	
	fp = open(sys.argv[1], 'r')

	mirna_ids = {}
	
	for line in fp:
		mirna_id = line.strip()

		if mirna_id not in mirna_ids:
			mirna_ids[mirna_id] = []

	#print (mirna_ids)

	new = []
	for dir in searchdirs:
		new.extend(list_mirnas(dir))


	# find matching ids

	for id in mirna_ids.keys():
		for case in new:
			if case.find(id+"_2")!=-1:
				mirna_ids[id].append(case)

	# print those with multiples

	for id in mirna_ids.keys():
		if len(mirna_ids[id]) > 0:

			for case in mirna_ids[id]:
				print("%s\t%s" % (id, case))

