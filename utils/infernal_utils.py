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

# --------------------------------------------------------------------------------------------------

import os
import sys
import re

START = 9
END = 10
STRAND = 11
EVAL = 2
BSCORE = 3
DBN_REGEX = {r"\(|\[|<|\{": '(', r"\)|\]|>|\}": ')', r"_|-|\.|,|~|:": '.'}

# --------------------------------------------------------------------------------------------------


def generate_bed_detail_file_with_ss(inf_output_file, dest_dir, ss_notation="wuss"):
    """
    Parses Infernal's detailed output and generates a bed file in detailed format
    with the last column containing the secondary structure string in the specified
    notation.

    inf_output_file: Infernal's output file (-o option)
    dest_dir: The path to the output directory
    ss_notation: A string indicating the the notation in which to output the
    secondary structure string (wuss or dbn)

    return: Void
    """

    scores = infernal_output_parser(inf_output_file, ss_notation=ss_notation)

    filename = os.path.basename(inf_output_file).partition('.')[0]
    fp_out = open(os.path.join(dest_dir, filename + '.bed'), 'w')

    for score in scores:
        # write output to file
        if score["strand"] == '+':
            fp_out.write(
                "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (score["rfamseq_acc"], score["start"], score["end"],
                                                      score["rna_type"], score["bit_score"], score["strand"],
                                                      score["e_value"], score["sec_struct"]))
        else:
            fp_out.write(
                "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (score["rfamseq_acc"], score["end"], score["start"],
                                                      score["rna_type"], score["bit_score"], score["strand"],
                                                      score["e_value"], score["sec_struct"]))
    fp_out.close()

# --------------------------------------------------------------------------------------------------


def infernal_output_parser(inf_output_file, ss_notation="wuss"):
    """
    Parses Infernal's detailed output file (-o) and returns a list of dictionaries
    for each hit found in the file

    inf_output_file: Infernal's output file (-o)
    dest_dir: The path to the output directory
    ss_notation: A string indicating the the notation in which to output the
    secondary structure string (wuss or dbn)

    return: A list of dictionaries
    """

    rna_type = ''
    ss_str_list = []

    scores = []

    fp_in = open(inf_output_file, 'r')

    line = fp_in.readline()
    rfam_acc = ''
    while line != '':

        # look for a new hit section and fetch rfam_acc
        if line.find("Accession:") != -1:
            rfam_acc = [x for x in line.strip().split(' ') if x != ''][1]

        # new hit
        if line[0:2] == ">>":
            line = line.strip()
            # get sequence accession
            seq_id = line[3:].strip()

            # get hit region
            index = 0
            while index < 3:
                line = fp_in.readline().strip()
                index += 1

            # bed file fields
            score_line = [x for x in line.strip().split(' ') if x != '']
            start = score_line[START]
            end = score_line[END]
            strand = score_line[11]
            e_value = score_line[2]
            bit_score = score_line[3]
            cm_start = score_line[6]
            cm_end = score_line[7]
            raw_trunc = score_line[14]
            truncated = ''

            if raw_trunc == 'no' or raw_trunc == '-':
                truncated = '0'
            elif raw_trunc.find('&') != -1:
                truncated = raw_trunc.replace('&', '').replace('\'', '')
            else:
                truncated = raw_trunc.replace('\'', '')

            rfamseq_acc = seq_id.split(' ')[0].split('|')[-1]

            # get secondary structure
            line = fp_in.readline()

            while line[0:2] != '>>' and line[0] != '#':

                if line.find("CS") != -1:
                    line = line.strip()
                    ss_line = [x for x in line.split(' ') if x != '' and x != 'CS']
                    ss_line = ss_line[0]
                    ss_str_list.append(ss_line)

                    if rna_type == '':  # 1st time we hit this
                        line = fp_in.readline().strip()  # read next line
                        rna_type = line.split(' ')[0].strip()

                line = fp_in.readline()  # read the rest

            sec_struct = ''.join(ss_str_list)

            if ss_notation.lower() == "dbn":
                sec_struct = convert_short_wuss_to_dbn(sec_struct)

            score_dict = {"rfam_acc": rfam_acc, "rfamseq_acc": rfamseq_acc, "start": start,
                          "end": end, "bit_score": bit_score, "e_value": e_value,
                          "strand": strand, "sec_struct": sec_struct, "rna_type": rna_type,
                          "truncated": truncated, "cm_start": cm_start, "cm_end": cm_end,
                          }

            scores.append(score_dict)

            ss_str_list = []
            rna_type = ''

        else:
            line = fp_in.readline()

    fp_in.close()

    return scores

# --------------------------------------------------------------------------------------------------


def infernal_to_rfam(inf_tblout_file, dest_dir, file_format='tsv'):
    """
    Parses Infernal's output file and exports results in Rfam's genome full region format
    (tsv option is used by default)

    inf_tblout_file: Infernal's output file in tabular format
    format: This is an option whether to output results in  tabular format or create a json file
    """

    in_file = open(inf_tblout_file, 'r')
    filename = os.path.basename(inf_tblout_file).partition('.')[0]

    out_file = None
    if file_format == "tsv":
        out_file = open(os.path.join(dest_dir, filename + ".tsv"), 'w')
    else:
        out_file = open(os.path.join(dest_dir, filename + ".json"), 'w')

    upid = filename
    is_significant = 1  # always the same before clan competition
    #seq_type = "null"
    seq_type = "full" # set all to full until we update seed sequences

    line = in_file.readline().strip()

    while line != '':
        if line[0] != '#':
            line = line.split(' ')
            line_cont = [x for x in line if x != '']

            genseq_acc = line_cont[2].split('|') # could also be position 0 for cmsearch
            genseq_acc = genseq_acc[len(genseq_acc) - 1]  # get last field
            rfam_acc = line_cont[1] # 1 or 3 for cmscan
            cm_start = line_cont[5]
            cm_end = line_cont[6]
            seq_start = line_cont[7]
            seq_end = line_cont[8]

            trunc = line_cont[10].replace('\'', '')
            if trunc == "no" or trunc == '-':
                trunc = str(0)

            elif trunc.find('&') != -1:
                trunc = trunc.replace('&', '')

            bit_score = line_cont[14]
            evalue = line_cont[15]

            out_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (rfam_acc, genseq_acc,
                                                                                 upid, seq_start,
                                                                                 seq_end, bit_score,
                                                                                 evalue, cm_start,
                                                                                 cm_end, trunc,
                                                                                 seq_type, is_significant))
            line = in_file.readline().strip()

        else:
            line = in_file.readline().strip()

    in_file.close()
    out_file.close()


# --------------------------------------------------------------------------------------------------


def tblout_to_full_region(tblout_file, dest_dir=None):
    """
    Parses Infernal's tblout file and generates a .txt file that is compatible with full_region
    table.

    tblout_file: A valid Infernal's output file in .tblout format
    dest_dir: The path to the output directory

    return: True if successful, False otherwise
    """

    tblout_fp = open(tblout_file, 'r')
    filename = os.path.split(tblout_file)[1].partition('.')[0]

    if dest_dir is None:
        dest_dir = os.path.split(tblout_file)[0]
        #filename = os.path.split(tblout_file)[1].partition('.')[0]

    full_region_fp = open(os.path.join(dest_dir, filename+'.txt'), 'w')

    # set to 1 on update and correct with clan competition
    is_significant = 1

    # set all to full until we update seed sequences
    seq_type = "full"

    line = tblout_fp.readline().strip()

    while line != '':
        # skip comment lines
        if str(line[0]) != '#':
            line = line.split(' ')
            line_cont = [x for x in line if x != '']
            if len(line_cont) > 1:
                seq_acc = ''
                if line[0].find('|') != -1:
                        seq_acc = line_cont[0].split('|')  # could also be position 0 for cmsearch
                        seq_acc = seq_acc[-1]  # get last field
                else:
                        seq_acc = line[0].strip()

                rfam_acc = line_cont[3]  # 1 or 3 for cmscan
                cm_start = line_cont[5]
                cm_end = line_cont[6]
                seq_start = line_cont[7]
                seq_end = line_cont[8]

                trunc = line_cont[10].replace('\'', '')
                if trunc == "no" or trunc == '-':
                        trunc = str(0)

                elif trunc.find('&') != -1:
                        trunc = trunc.replace('&', '')

                bit_score = line_cont[14]
                evalue = line_cont[15]

                full_region_fp.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (rfam_acc,
                                                                                   seq_acc,
                                                                                   seq_start,
                                                                                   seq_end,
                                                                                   bit_score,
                                                                                   evalue,
                                                                                   cm_start,
                                                                                   cm_end,
                                                                                   trunc,
                                                                                   seq_type,
                                                                                   str(is_significant)))
            line = tblout_fp.readline().strip()

        else:
            line = tblout_fp.readline().strip()

    tblout_fp.close()
    full_region_fp.close()


# --------------------------------------------------------------------------------------------------


def infernal_to_full_region(inf_output_file, dest_dir, filename=None):
    """
    Parses Inferna's detailed output (-o option) and generates a file in tabular format, which is
    compatible with the full_region table

    inf_output_file: Infernal's output file (-o option)
    dest_dir: The path to the output directory
    filename: A filename for the output

    returns: Void
    """

    scores = infernal_output_parser(inf_output_file, ss_notation="wuss")

    if dest_dir is None:
        dest_dir = os.path.split(inf_output_file)[0]

    if filename is None:
        filename = os.path.basename(inf_output_file).partition('.')[0]

    fp_out = open(os.path.join(dest_dir, filename + '.txt'), 'w')

    for score in scores:
        fp_out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (score["rfam_acc"],
                                                                       score["rfamseq_acc"],
                                                                       score["start"],
                                                                       score["end"],
                                                                       score["bit_score"],
                                                                       score["e_value"],
                                                                       score["cm_start"],
                                                                       score["cm_end"],
                                                                       score["truncated"],
                                                                       "full", '1'))
    fp_out.close()


# --------------------------------------------------------------------------------------------------


def convert_short_wuss_to_dbn(ss_string):
    """
    Converts RNA structure string from shorthand WUSS notation to dot-bracket notation

    ss_string: Secondary structure string
    """

    for regex in DBN_REGEX.keys():
        ss_string = re.sub(regex, DBN_REGEX[regex], ss_string)

    return ss_string


# --------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    pass



