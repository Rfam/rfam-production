import os
import sys
import argparse

from subprocess import Popen, PIPE

# ---------------------------------------------------------------------------------------------


def validate_xml_dump(xml_file):
    """
    Validates Rfam XML dumps using xmllint

    xml_file: An XML dump used to index Rfam release data

    return: True upon success, False for failure
    """
    process = Popen(["xmllint", "--schema", "http://www.ebi.ac.uk/ebisearch/XML4dbDumps.xsd", "--noout", xml_file],
                    stdin=PIPE, stdout=PIPE, stderr=PIPE)
    output = process.communicate()[1]

    # check if output contains the keyword validates
    if output.find("validates") == -1:
        return False

    return True


# ---------------------------------------------------------------------------------------------


def parse_arguments():
    """
    Simple argument parsing using python's argparse

    return: Python's argparse parser object
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("--input", help="Single XML file or directory", action="store")
    parser.add_argument("--log", help="Generate a log file listing all XML files failining validation",
                        action="store_true")

    return parser

# ---------------------------------------------------------------------------------------------


if __name__ == '__main__':

    parser_obj = parse_arguments()
    args = parser_obj.parse_args()

    if os.path.isfile(args.input):
        if validate_xml_dump(sys.argv[1]) is False:
            print (args.input)

    elif os.path.isdir(args.input):
        fp = None

        if args.log:
            fp = open(os.path.join(args.input, "error.log"), 'w')

        xml_files = [x for x in os.listdir(args.input)]
        for xml_file in xml_files:
            if validate_xml_dump(os.path.join(args.input, xml_file)) is False:
                if args.log:
                    fp.write(xml_file + '\n')
                else:
                    print (xml_file)

        if args.log:
            fp.close()
