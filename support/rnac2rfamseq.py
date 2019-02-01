import os
import sys
import subprocess

from config import rfam_config as rc

# ----------------------------------------------------------------------------------


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
        tax_id = rnac_line_contents[3]

        # append taxid to URS to create RNAcentral accession
        rnac_acc = rnac_line_contents[0] + '_' + tax_id
        # all URS entries must be unique, count if any skipped
        if rnac_acc not in rnac_metadata_dict:
            rnac_metadata_dict[rnac_acc] = {}
            rnac_metadata_dict[rnac_acc]["source"] = rnac_line_contents[1]
            temp = rnac_line_contents[2].partition('/')

            previous_acc = temp[0]
            """ Revise this section
            if temp[2] != '':
                region_coords = temp[2].split('-')
                seq_start = region_coords[0]
                seq_end = region_coords[1]
            else:
                seq_start = 1
                seq_end = seq_length
            """
            rnac_metadata_dict[rnac_acc]["previous_acc"] = previous_acc
            rnac_metadata_dict[rnac_acc]["tax_id"] = tax_id
            rnac_metadata_dict[rnac_acc]["mol_type"] = rnac_line_contents[4]

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

# ----------------------------------------------------------------------------------


def generate_sequence_stats(fasta_file):
    """
    Uses esl-sfetch to generate stats for every sequence in the provided fasta file
    and loads that information in a dictionary

    fasta_file: A sequence file in fasta format

    return: A dictionary with sequence statistics
    """

    sequence_stats = {}
    # create argument list required by subprocess PIPE
    args = [rc.ESL_SEQSTAT, "-a", "--dna", fasta_file]

    # open a process pipe and run esl-seqstat. Stores the result in process
    process = subprocess.Popen(args, stdout=subprocess.PIPE)

    # fetch esl-seqstat results
    result = process.communicate()

    # read in and clean-up sequence stats
    seq_stats = [x[2:] for x in result[0].strip().split('\n') if x[0] == '=']

    # process every sequence header to generate a rfamseq entry
    for sequence in seq_stats:
        seq_elements = [x.strip() for x in sequence.split(' ') if x != '']

        rfamseq_acc = seq_elements[0]

        if rfamseq_acc not in sequence_stats:
            sequence_stats[rfamseq_acc] = {}
            sequence_stats[rfamseq_acc]["length"] = seq_elements[1]
            sequence_stats[rfamseq_acc]["desc"] = ' '.join(seq_elements[2:])

    return sequence_stats

# ----------------------------------------------------------------------------------


def generate_rfamseq_tsv_dump(fasta_file, rnac_metadata_file, to_file = False):
    """
    Parses the sequence file and the rnacentral metadata file and generates an
    rfamseq dump for the sequences to be imported to rfam_live

    fasta_file: The path to the input fasta file as generated from RNAcentral
    rnac_metadata_file: The path to the RNAcentral metadata file in tabular
    format

    return: void
    """

    rnac_metadata = load_rnacentral_metadata_to_dict(rnac_metadata_file)
    sequence_stats = generate_sequence_stats(fasta_file)

    for urs_id in rnac_metadata.keys():
        rfamseq_acc = urs_id
        accession = urs_id
        version = '0'
        taxid = rnac_metadata[urs_id]["tax_id"]
        mol_type = rnac_metadata[urs_id]["mol_type"]
        length = sequence_stats[urs_id]["length"]
        description = sequence_stats[urs_id]["desc"]
        previous_acc = rnac_metadata[urs_id]["previous_acc"]
        source = rnac_metadata[urs_id]["source"]


        rfamseq_entry = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (rfamseq_acc, accession,
                                                          version, taxid,
                                                          mol_type, length,
                                                          description, previous_acc,
                                                          source)

        # TO DO provide the option to generate a file instead?

        print rfamseq_entry

# ----------------------------------------------------------------------------------

if __name__ == "__main__":


    fasta_file = sys.argv[1]
    rnac_metadata_file = sys.argv[2]

    generate_rfamseq_tsv_dump(fasta_file, rnac_metadata_file, to_file=False)
