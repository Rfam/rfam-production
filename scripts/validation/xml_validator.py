import os
import sys

from subprocess import Popen, PIPE

# ---------------------------------------------------------------------------------------------


def validate_xml_dump(xml_file):
	
	process = Popen(["xmllint", "--schema", "http://www.ebi.ac.uk/ebisearch/XML4dbDumps.xsd", "--noout", xml_file], stdin=PIPE, stdout=PIPE, stderr=PIPE)
        output = process.communicate()[1]

	# check if output contains the keyword validates
    	if output.find("validates") == -1:
        	return False

    	return True

# ---------------------------------------------------------------------------------------------

if __name__=='__main__':

	if os.path.isfile(sys.argv[1]):
		if validate_xml_dump(sys.argv[1]) is False:
			print (sys.argv[1])

	elif os.path.isdir(sys.argv[1]):
		xml_files = [x for x in os.listdir(sys.argv[1])]
		for xml_file in xml_files:
			if validate_xml_dump(os.path.join(sys.argv[1], xml_file)) is False:
				print (xml_file)
