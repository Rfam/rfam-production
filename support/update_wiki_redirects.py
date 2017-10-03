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

# ---------------------------------IMPORTS-------------------------------------

import os
import sys
import subprocess
import logging


# -----------------------------------------------------------------------------

def wiki_redirects_parser(redirects):
    """
    Parses wiki redirects and reports rfam accessions and changes to WK tags
    Returns a dictionary with changes per family accession

    redirects: redirects.txt file from wiki
    """

    wk_edits = {}
    fp = open(redirects, 'r')

    rfam_lines = [x.strip() for x in fp if x.find('Rfam') != -1]
    fp.close()

    c = '"'
    for line in rfam_lines:
        rfam_acc = line[len(line) - 7:len(line)]
        positions = [pos for pos, char in enumerate(line) if char == c]
        old_str = line[positions[0] + 1:positions[1]]
        new_str = line[positions[2] + 1:positions[3]]
        wk_edits[rfam_acc] = (old_str, new_str)

    # return a dictionary with the new redirects per family
    return wk_edits


# -----------------------------------------------------------------------------

def update_desc_file(desc_file, tag, updates):
    """
    Function to update changes in desc files automatically

    desc_file: Path to a desc file
    tag: A valid DESC file tag to be modified e.g. 'WK'
    updates: A tuple with the desc changes (old, new)
    """
    family_dir = os.path.split(desc_file)[0]
    new_desc_path = os.path.join(family_dir, "DESC_NEW")

    fp = open(desc_file, 'r')
    new_desc = open(new_desc_path, 'w')
    desc_lines = fp.readlines()
    fp.close()

    for line in desc_lines:
        if line.find(tag) != -1:
            line = line.replace('_', ' ')
            line = line.replace(updates[0], updates[1].replace(' ', '_'))
        new_desc.write(line)

    new_desc.close()

    # remove old desc and rename new one to DESC
    os.remove(desc_file)
    os.rename(new_desc_path, os.path.join(family_dir, "DESC"))


# -----------------------------------------------------------------------------

def checkout_family_from_svn(rfam_acc, dest_dir):
    """
    Checks out a family from the svn by calling rfco

    rfam_acc: A valid Rfam family accession
    dest_dir: Destination directory where to check out family
    """

    os.chdir(dest_dir)
    cmd = "rfco %s" % (rfam_acc)
    subprocess.call(cmd, shell=True)

    if not os.path.exists(os.path.join(dest_dir, rfam_acc)):
        return -1

    return 0


# -----------------------------------------------------------------------------

def check_family_into_svn(dest_dir, onlydesc=None):
    """
    Check family back to SVN repo using rfci

    :param dest_dir:
    :return: void
    """

    os.chdir(dest_dir)

    # list all family directories
    family_dirs = [x for x in os.listdir(dest_dir) if x.find('RF') != -1]

    for family in family_dirs:
        if onlydesc is None:
            cmd = "rfci -m /'Updated WK in DESC/' %s" % family
        else:
            cmd = "rfci -onlydesc -m /'Updated WK in DESC/' %s" % family

        subprocess.call(cmd, shell=True)


# -----------------------------------------------------------------------------

def create_wiki_markdown_links(dest_dir):
    """
    Create a list of Rfam accessions and links in markdown

    dest_dir: The check out directory
    """

    fp_out = open(os.path.join(dest_dir, "family_links.md"), 'w')
    families = [x for x in os.listdir(dest_dir) if os.path.isdir(os.path.join(dest_dir, x))]

    for rfam_acc in families:
        fam_dir = os.path.join(dest_dir, rfam_acc)
        desc_fp = open(os.path.join(fam_dir, "DESC"), 'r')
        description = ''
        wk_tag = ''

        # get tags
        for line in desc_fp:
            if line[0:2] == 'DE':
                description = line[2:].strip()
            elif line[0:2] == 'WK':
                wk_tag = line[2:].replace(' ', '')

        desc_fp.close()

        md_str = ''
        md_str = '[' + rfam_acc + '-' + description + ']'
        md_str = md_str + "(https://en.wikipedia.org/wiki/" + wk_tag + ')'
        fp_out.write(md_str + '\n')

    fp_out.close()


# -----------------------------------------------------------------------------

def create_rfam_markdown_links(dest_dir):
    """
    Create a list of Rfam accessions and links in markdown

    dest_dir: The check out directory
    """

    fp_out = open(os.path.join(dest_dir, "family_links.md"), 'w')
    families = [x for x in os.listdir(dest_dir) if os.path.isdir(os.path.join(dest_dir, x))]

    for rfam_acc in families:
        fam_dir = os.path.join(dest_dir, rfam_acc)
        desc_fp = open(os.path.join(fam_dir, "DESC"), 'r')
        wk_tag = ''

        # get tags
        for line in desc_fp:
            if line[0:2] == 'WK':
                wk_tag = line[2:].replace(' ', '')
        desc_fp.close()

        fp_out.write('[' + rfam_acc + '-' + wk_tag + ']' + "(http://rfam.xfam.org/family/" + rfam_acc + ")\n\n")

    fp_out.close()


# -----------------------------------------------------------------------------

def main(redirects, dest_dir):
    """
    This function puts all the pieces together parameters are provided through
    command line

    redirects: Wiki redirects output file
    dest_dir: Family check out directory
    """

    # create a log file
    logging.basicConfig(filename=os.path.join(dest_dir, 'wk_desc_updates.log'), level=logging.ERROR)

    # parse wiki redirects file and changes per rfam_acc in a dictionary
    wk_edits = wiki_redirects_parser(redirects)

    for rfam_acc in wk_edits.keys():
        # check out family from the svn repo
        status = checkout_family_from_svn(rfam_acc, dest_dir)
        # rfco success
        if status == 0:
            family_dir = os.path.join(dest_dir, rfam_acc)
            update_desc_file(os.path.join(family_dir, "DESC"), "WK", wk_edits[rfam_acc])
        else:
            logging.error("Failed to checkout family %s" % rfam_acc)

    # generate markdown for inspection
    create_rfam_markdown_links(dest_dir)


# --------------------------------------------------------------------------------------------------

def commit_family_to_svn(dest_dir):
    """
    Commits a list of families to the svn using rfci -onlydesc

    dest_dir: A destination directory with all
    """

    os.chdir(dest_dir)

    # create a log file
    logging.basicConfig(filename=os.path.join(dest_dir, 'auto_rfci_errors.log'), level=logging.DEBUG)

    family_dirs = [x for x in os.listdir(dest_dir) if os.path.isdir(os.path.join(dest_dir, x))]

    rfci_cmd = "rfci -onlydesc -m \'Wiki Updates\' %s"

    for family_dir in family_dirs:
        cmd = rfci_cmd % family_dir

        try:
            subprocess.call(cmd, shell=True)

        except:
            logging.exception("Failed to commit family %s" % family_dir)


# --------------------------------------------------------------------------------------------------

def usage():
    # CODE
    pass


# --------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    redirects_file = sys.argv[1]
    dest_dir = sys.argv[2]

    main(redirects_file, dest_dir)

    # need to wrap this up in an option to just commit or develop an autocommit script
    # commit_family_to_svn(dest_dir)
