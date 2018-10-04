import os
import sys

"""
Script to convert uniprot proteome files to rfam pipeline files such genome_downloader,
genome_scanner etc.
"""

# ----------------------------------------------------------------------------------


def convert_proteome_list_to_upid_gca_file(input_file, filename=None, dest_dir=None):
    """
    Converts a proteome list file in tabular format to the upid_gca file format
    used with the genome_downloader and genome_search pipelines

    input_file: A proteome-all.tab file downloaded from Uniprot's proteomes website
    filename: A name for the new file to be generated. Defaults to input_file filename
    if None
    dest_dir: A path to an output directory. If None, uses the location of the input
    file
    """

    input_fp = open(input_file, 'r')

    if filename is None:
        filename = os.path.basename(input_file).partition('.')[0]

    if dest_dir is None:
        dest_dir = os.path.split(input_file)[0]

    output_fp = open(os.path.join(dest_dir, filename + '.tsv'), 'w')

    for proteome in input_fp:
        fields = proteome.strip().split('\t')

        # mark any erroneous cases with  missing domains
        if len(fields) < 4 and fields[-1].find(',') == -1:
            domain = "N/A"
        else:
            domain = fields[3].split(',')[0].strip().lower()

        output_fp.write(fields[0] + '\t' + fields[2] + '\t' + domain + '\n')

    output_fp.close()


# ----------------------------------------------------------------------------------


def dump_new_upids_to_file(upid_gca_old, upid_gca_new, dest_dir=None):
    """
    Compares two upid_gca files from different Uniprot proteome releases
    and dumps the new upids in a new upid_gca file all resulting from the
    newer version of the proteome file (upid_gca_new). To be used for
    downloading a subset of genomes.

    upid_gca_out: A valid upid_gca file in .tsv format of an older reference
    proteome release
    upid_gca_new: A valid upid_gca file in .tsv format of a newer reference
        proteome release
    dest_dir: The path to a destination directory. Uses the directory of
    upid_gca_new if dest_dir is None

    returns: void
    """

    # open and load the accessions to a dictionary
    old_upids = {}
    upid_gca_old_fp = open(upid_gca_old, 'r')
    upid_gca_new_fp = open(upid_gca_new, 'r')

    # load old accessions
    for line in upid_gca_old_fp:
        fields = line.strip().split('\t')
        if fields[0] not in old_upids:
            old_upids[fields[0]] = fields[1:]

    upid_gca_old_fp.close()

    # loop over new upids
    upid_gca_new_fp = open(upid_gca_new, 'r')

    # if no destination directory provided, use directory of upid_gca_new
    if dest_dir is None:
        dest_dir = os.path.split(upid_gca_new)[0]

    upid_gca_out_fp = open(os.path.join(dest_dir, "new_upids.tsv"), 'w')

    for line in upid_gca_new_fp:
        fields = line.strip().split('\t')
        # check if upid in the old release
        if fields[0] not in old_upids:
            upid_gca_out_fp.write(('\t').join(fields) + '\n')

    upid_gca_out_fp.close()


# ----------------------------------------------------------------------------------

if __name__ == '__main__':

    if "--new-upids" in sys.argv:
        upid_gca_old = sys.argv[1]
        upid_gca_new = sys.argv[2]

        if len(sys.argv) == 5:
            dest_dir = sys.argv[3]
            dump_new_upids_to_file(upid_gca_old, upid_gca_new, dest_dir)

        elif len(sys.argv) == 4:
            dump_new_upids_to_file(upid_gca_old, upid_gca_new, None)

    elif "--upid-gca-file" in sys.argv:
        if len(sys.argv) == 5:
            input_file = sys.argv[1]  # proteome.tab file dowloaded from Uniprot
            filename = sys.argv[2]
            dest_dir = sys.argv[3]
            convert_proteome_list_to_upid_gca_file(input_file, filename, dest_dir)

        elif len(sys.argv) == 4:
            input_file = sys.argv[1]  # proteome.tab file dowloaded from Uniprot
            filename = sys.argv[2]
            convert_proteome_list_to_upid_gca_file(input_file, filename, None)

        elif len(sys.argv) == 3:
            input_file = sys.argv[1]  # proteome.tab file dowloaded from Uniprot
            convert_proteome_list_to_upid_gca_file(input_file, None, None)

    else:
        # need to develop a help function and use arg parse?
        sys.exit("Wrong parameters! Please try again")
