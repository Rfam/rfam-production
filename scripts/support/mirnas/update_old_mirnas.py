import os
import sys
import json
import shutil
import subprocess

searchdirs = ["/hps/nobackup/production/xfam/rfam/RELEASES/14.3/miRNA_relabelled/batch1_chunk1_searches",
              "/hps/nobackup/production/xfam/rfam/RELEASES/14.3/miRNA_relabelled/batch1_chunk2_searches",
              "/hps/nobackup/production/xfam/rfam/RELEASES/14.3/miRNA_relabelled/batch2/searches"]

# ------------------------------------------------------------------------


def checkout_family(rfam_acc, dest_dir):

	os.chdir(dest_dir)
	cmd = "rfco.pl %s" % rfam_acc

	subprocess.call(cmd, shell=True)

	if os.path.exists(os.path.join(dest_dir, rfam_acc)):
		return True

	return False


# ------------------------------------------------------------------------


if __name__ == "__main__":

	dest_dir = "/hps/nobackup/production/xfam/rfam/RELEASES/14.3/miRNA_relabelled/update_old_rfam_mirnas"

	mirna_families = sys.argv[1]

	fp = open(mirna_families, 'r')

	mirna_dict = json.load(fp)

	fp.close()


	mirnas_to_make = {}

	for mirna in mirna_dict:
		mirna_dir = mirna
		if mirna.find("relabelled")==-1:
			mirna_dir = mirna+'_relabelled'
		
		if mirna_dict[mirna]["threshold"] != '' and mirna_dict[mirna]["rfam_match"] == 100:
			rfam_acc = mirna_dict[mirna]["rfam_acc"]
			mirnas_to_make[rfam_acc] = ""
			#check = checkout_family(rfam_acc, dest_dir)
			checkout_dir = os.path.join(dest_dir, rfam_acc)
			check = True
			if check is True:
				updated_mirna_loc = ""
				for searchdir in searchdirs:
					if os.path.exists(os.path.join(searchdir, mirna_dir)):
						updated_mirna_loc = os.path.join(searchdir, mirna_dir)
						updated_seed = os.path.join(updated_mirna_loc, "SEED")
						os.rename(os.path.join(checkout_dir, "SEED"), os.path.join(checkout_dir, "SEED_old"))
						shutil.copyfile(updated_seed, os.path.join(checkout_dir, "SEED"))					
					else:
						continue

	fp = open("old_mirnas_to_make.json", 'w')
	json.dump(mirnas_to_make, fp)
	fp.close()
