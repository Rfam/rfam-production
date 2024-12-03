import os
import sys

if __name__=="__main__":

	locations_file = sys.argv[1]

	fp = open(locations_file, 'r')
	fam_dirs = [x.strip() for x in fp]
	fp.close()


	for dir_loc in fam_dirs:
		err_file_loc = os.path.join(dir_loc, "auto_rfmake.err")
		if os.path.exists(err_file_loc):
			os.remove(err_file_loc)
		out_file_loc = os.path.join(dir_loc, "auto_rfmake.out")
		if os.path.exists(out_file_loc):
			os.remove(out_file_loc)


	print "Done"
