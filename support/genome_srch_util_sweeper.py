import os
import sys
import shutil

if __name__=='__main__':
	project_dir = sys.argv[1]
	tool = sys.argv[2]

	subdirs = [x for x in os.listdir(project_dir) if os.path.isdir(os.path.join(project_dir, x))]
	
	for subdir in subdirs:
		subdir_loc = os.path.join(project_dir, subdir)
		updirs = os.listdir(subdir_loc)
		
		for updir in updirs:
			updir_loc = os.path.join(subdir_loc, updir)
			
			if tool.lower()=='merge':

				merge_err = os.path.join(updir_loc, "merge.err")
				if os.path.exists(merge_err):
					os.remove(merge_err)
			
				merge_out = os.path.join(updir_loc, "merge.out")
				if os.path.exists(merge_out):
					os.remove(merge_out)

			elif tool.lower()=='split':
				chunk_dir = os.path.join(updir_loc, "search_chunks")
				if os.path.exists(chunk_dir):
					shutil.rmtree(chunk_dir)

				split_err = os.path.join(updir_loc, "split.err")
				if os.path.exists(split_err):
                                        os.remove(split_err)
				
				split_out = os.path.join(updir_loc, "split.out")
				if os.path.exists(split_out):
                                        os.remove(split_out)

			elif tool.lower()=='scan':
				output_dir = chunk_dir = os.path.join(updir_loc, "search_output")
				if os.path.exists(output_dir):
					shutil.rmtree(output_dir)
			
