import os
import sys

from config import rfam_config as rc

# ----------------------------------------------------------------------

def load_rnacentral_metadata_to_dict(rnac_metadata_file, to_file=False, destdir=None):
    """
    Loads RNAcentral metadata to a dictionary. The expected format is:
    URS_acc\tsource\tpre_acc/start-end\tncbi_id\tmol_type

    rnac_metadata_file: A tab delimited file with Rfam compatible metadata
    extracted from RNAcentral
    to_file: If True, it dumps the output dictionary in a json file
    destdir: The path to a destination directory where any output will be generated

    return: A dictionary with all metadata values loaded from rnacentral tsv file
    """

    # open a new file handle
    fp = open(rnac_metadata_file, 'r')

    rnac_metadata_dict = {}
    # loop over all lines in the rnacentral tsv file
    for rnac_line in fp:
        rnac_line_contents = rnac_line.strip().split('\t')
        # all URS entries must be unique, count if any skipped
        if rnac_line_contents[0] not in rnac_metadata_dict:
            rnac_metadata_dict[rnac_line_contents[0]] = {}
            rnac_metadata_dict[rnac_line_contents[0]]["source"] = rnac_line_contents[1]
            seq_length = rnac_line_contents[3]
            temp = rnac_line_contents[2].partition('/')

            previous_acc = temp[0]
            if temp[2] != '':
                region_coords = temp[2].split('-')
                seq_start = region_coords[0]
                seq_end = region_coords[1]
            else:
                seq_start = 1
                seq_end = seq_length

            rnac_metadata_dict[rnac_line_contents[0]]["previous_acc"] = previous_acc
            rnac_metadata_dict[rnac_line_contents[0]]["seq_start"] = seq_start
            rnac_metadata_dict[rnac_line_contents[0]]["seq_end"] = seq_end
            rnac_metadata_dict[rnac_line_contents[0]]["seq_length"] = seq_length
            rnac_metadata_dict[rnac_line_contents[0]]["mol_type"] = rnac_line_contents[4]

    # close file handle
    fp.close()

    # create a json dump
    if to_file is True:
        import json
        filename = os.path.basename(rnac_metadata_file).partition('.')[0] + '.json'

        if destdir is None:
            destdir = os.path.split(rnac_metadata_file)[0]

        json_dump_file = os.path.join(os.path.join(destdir, filename))

        fp = open(json_dump_file, 'w')
        json.dump(rnac_metadata_dict, fp)

    return rnac_metadata_dict

# ----------------------------------------------------------------------

if __name__ == "__main__":

    pass


