"""
Copyright [2009-2017] EMBL-European Bioinformatics Institute
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

# --------------------------------------------------------------------------------

import os
import sys
import subprocess

# --------------------------------------------------------------------------------


def parse_go_term_file(go_term_file):
    """
    Parses a tabular file of valid GO terms and returns a dictionary with all
    available GO terms per family. An example of the file can be found in this
    url(https://goo.gl/5PN5LT)

    go_term_file: The path to a go term file. This should be a tab delimited
    file where multiple terms followed by their description, are

    return: void
    """

    # need to add a column for the function DELETE, ADD, REPLACE

    go_terms_dict = {}
    go_terms_fp = open(go_term_file, 'r')

    for go_line in go_terms_fp:
        go_line = go_line.strip().split('\t')

        go_terms = go_line[1:]

        rfam_acc = go_line[0].strip()

        if rfam_acc not in go_terms_dict:
            go_terms_dict[rfam_acc] = []

        # check if the go term pairs are correct. Each GO term id must be followed
        # by a description
        if len(go_terms) % 2 != 0:
            print "Erroneous go term for %s" % rfam_acc
            continue

        idx = 0
        # work on go term pairs for a specific family
        while idx < len(go_terms):
            # organize go term pairs in tuples
            go_terms_dict[rfam_acc].append((go_terms[idx].strip(),
                                            go_terms[idx+1].strip()))
            idx += 2

    return go_terms_dict

# --------------------------------------------------------------------------------


def parse_go_term_validation_file(go_validation_file):
    """


    param go_validation_file:
    returns: A dict
    """

    go_updates = {}

    val_file_fp = open(go_validation_file, 'r')

    for err_line in val_file_fp:
        err_line = err_line.strip().split(':')[1:]
        rfam_acc = err_line[0]

        # check first if a new family line
        if rfam_acc[0:2] == "RF":
            # create a list for the new family
            if rfam_acc not in go_updates.keys():
                go_updates[rfam_acc] = {}

            error_text = err_line[2]

            # check for error type
            if error_text.find("OBSOLETE") != -1 or error_text.find("UNKNOWN") != -1:
                go_term_id = error_text.split(' ')[0].strip()
                go_updates[rfam_acc].update({go_term_id: {"FUNCTION": "DELETE"}})

            elif error_text.find("SECONDARY") != -1:
                new_go = err_line[len(err_line)-1]
                # replace or do nothing..
                #go_updates[rfam_acc].append("SECONDARY")

            elif error_text.find("incorrect") != -1:
                go_term_id = error_text.split(' ')[0].strip()
                replacement_text = error_text.partition("-")[2].strip().partition("should be")
                old_text = replacement_text[0].replace('\"', '').strip()
                new_text = replacement_text[2].replace('\"', '').strip()
                go_updates[rfam_acc].update({go_term_id: {"OLD": old_text,
                                             "NEW": new_text, "FUNCTION": "REPLACE"}})

    return go_updates

# --------------------------------------------------------------------------------


def update_desc_go_terms(go_term_list, desc_file_path):
    """
    Modifies a valid family DESC file, by adding GO terms listed in go_term_list

    go_term_list: A list of GO term/description tuples to include in the DESC file
    desc_file_path: The path to a valid Rfam DESC file

    return: void
    """

    desc_fp = open(desc_file_path, 'r')
    desc_lines = desc_fp.readlines()
    desc_fp.close()

    new_desc_fp = open(os.path.join(os.path.split(desc_file_path)[0], "DESC_new"), 'w')

    idx = 0
    # loop once until we reach DR lines
    while idx < len(desc_lines):

        # break the loop if DR line found
        if desc_lines[idx][0:2] == 'DR':
            break

        new_desc_fp.write(desc_lines[idx])
        idx += 1

    # loop over all DR lines
    while idx < len(desc_lines):
        # break the loop if DR line found
        if desc_lines[idx][0:2] != 'DR':
            break

        # skip all GO terms
        elif desc_lines[idx].find('GO') == -1:
            new_desc_fp.write(desc_lines[idx])

        idx += 1

    # add new GO lines (e.g. GO; 0006396; RNA processing;)
    for go_pair in go_term_list:
        new_desc_fp.write('DR   GO; ' + go_pair[0].replace('GO:', '') + '; ' + go_pair[1] + ';\n')

    # write the remaining DESC lines
    while idx < len(desc_lines):
        new_desc_fp.write(desc_lines[idx])

        idx += 1

    new_desc_fp.close()

    # delete old DESC
    os.remove(desc_file_path)

    # rename DESC_new to DESC
    os.rename(os.path.join(os.path.split(desc_file_path)[0], "DESC_new"),
              os.path.join(os.path.split(desc_file_path)[0], "DESC"))

# --------------------------------------------------------------------------------


def modify_desc_go_terms(go_term_dict, desc_file_path):
    """
    List of GO term modifications for a specific family parse_go_term_validation_file
    The function reads the

    go_term_list: Family GO term list as generated by
    desc_file_path:
    :return:
    """

    # need to do the check for an empty list outside this function

    # get destination for new desc

    # read desc lines and delete DESC file
    desc_fp = open(desc_file_path, 'r')
    desc_lines = desc_fp.readlines()
    desc_fp.close()
    os.remove(desc_file_path)

    new_desc_fp = open(desc_file_path, 'w')

    for desc_line in desc_lines:
        if desc_line[0:2] != "DR":
            new_desc_fp.write(desc_line)
        else:
            # get go_term
            desc_go_term = desc_line.split(';')[1].strip()
            # if GO term is not in dictionary, write it to file as we are keeping it

            if desc_go_term not in go_term_dict.keys():
                new_desc_fp.write(desc_line)
            else:
                # check if empty and skip
                if bool(go_term_dict) is False:
                    continue

                # if go term is marked for deletion move to next iteration
                elif go_term_dict[desc_go_term]["FUNCTION"] == "DELETE":
                    continue

                elif go_term_dict[desc_go_term]["FUNCTION"] == "REPLACE":
                    new_go_line = desc_line.replace(go_term_dict[desc_go_term]["OLD"],
                                                    go_term_dict[desc_go_term]["NEW"])
                    new_desc_fp.write(new_go_line)

                elif go_term_dict[desc_go_term]["FUNCTION"] == "ADD":
                    # to be implemented
                    pass

    new_desc_fp.close()

# --------------------------------------------------------------------------------


def main(go_term_file, option ,checkout_dir):
    """
    This function will check out a family from the SVN repository, modify the
    desc file to add the GO terms and then check family back to the SVN repository

    go_term_file: The path to a go term file. This should be a tab delimited
    file where multiple terms followed by their description, are
    checkout_dir: A destination directory where to checkout all families

    return: void
    """

    checkout_cmd = "rfco.pl %s"
    checkin_cmd = "rfci.pl -onlydesc -m \'GO terms added\' %s"

    if not os.path.exists(checkout_dir):
        os.mkdir(checkout_dir)

    os.chdir(checkout_dir)

    # parse go term file
    go_term_dict = None
    if option == "--curation":
        go_term_dict = parse_go_term_file(go_term_file)

    elif option == "--validation":
        go_term_dict = parse_go_term_validation_file(go_term_file)

    for rfam_acc in go_term_dict.keys():
        fam_go_terms = go_term_dict[rfam_acc]

        # make upper case
        rfam_acc = rfam_acc.upper()

        # checkout family if go term dictionary is not empty
        if bool(fam_go_terms) is not False:
            subprocess.call(checkout_cmd % rfam_acc.upper(), shell=True)

            family_dir = os.path.join(checkout_dir, rfam_acc.upper())

            # add GO terms to DESC
            if option == "--curation":
                update_desc_go_terms(fam_go_terms, os.path.join(family_dir, 'DESC'))

            elif option == "--validation":
                modify_desc_go_terms(fam_go_terms, os.path.join(family_dir, 'DESC'))

            # checkin family to SVN
            subprocess.call(checkin_cmd % rfam_acc.upper(), shell=True)


# --------------------------------------------------------------------------------

if __name__ == '__main__':

    go_term_file = sys.argv[1]
    checkout_dir = sys.argv[2]
    option = sys.argv[3]

    main(go_term_file, option, checkout_dir)
