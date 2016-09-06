"""
Copyright [2009-2016] EMBL-European Bioinformatics Institute
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
import re

START = 9
END = 10
STRAND = 11
EVAL = 2
BSCORE = 3
DBN_LEFT = '('
DBN_REGEX = {r"\(|\[|<|\{": '(', r"\)|\]|>|\}": ')', r"_|-|\.|,|~|:": '.'}


# --------------------------------------------------------------------------------------------------

def infernal_output_parser(inf_output_file, dest_dir, ss_notation="wuss"):
    """
    Exports secondary structure from infernal's output file.

    inf_output_file: Infernal's output file
    """

    rna_type = ''
    ss_str_list = []

    fp = open(inf_output_file, 'r')

    out_file_name = os.path.basename(inf_output_file)

    if out_file_name.find('.') != -1:
        out_file_name = out_file_name.partition('.')[0]

    fp_out = open(os.path.join(dest_dir, out_file_name + ".bed"), 'w')

    line = fp.readline()

    while line != '':
        # new hit
        if line[0:2] == ">>":
            line = line.strip()
            # get sequence accession
            seq_id = line[3:].strip()

            # get hit region
            index = 0
            while index < 3:
                line = fp.readline().strip()
                index += 1

            # bed file fields
            score_line = filter(lambda x: x != '', line.strip().split(' '))

            start = score_line[START]
            end = score_line[END]
            strand = score_line[11]
            e_value = score_line[2]
            bit_score = score_line[3]

            # get secondary structure
            line = fp.readline()

            while line[0:2] != '>>' and line[0] != '#':

                if line.find("CS") != -1:
                    line = line.strip()
                    ss_line = [x for x in line.split(' ') if x != '' and x != 'CS']
                    ss_line = ss_line[0]
                    ss_str_list.append(ss_line)

                    if rna_type == '':  # 1st time we hit this
                        line = fp.readline().strip()  # read next line
                        rna_type = line.split(' ')[0].strip()

                line = fp.readline()  # read the rest

            sec_struct = ''.join(ss_str_list)

            if ss_notation.lower() == "dbn":
                sec_struct = convert_short_wuss_to_dbn(sec_struct)

            # write output to file
            fp_out.write(
                "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (seq_id, start, end, strand,
                                                      rna_type, e_value, bit_score, sec_struct))

            ss_str_list = []

        else:
            line = fp.readline()

    fp.close()
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

