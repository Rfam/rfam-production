import os
import sys
import json
import subprocess
from config import rfam_config as gf

# -----------------------------------------------------------------------------


def generate_rfamseq_metadata_from_fasta(fasta_file, taxid, source, filename=None, to_file=True):
    """
    Parses a fasta file and generates rfamseq like matadata using esl-seqstat. The output is in
    rfamseq table format.

    fasta_file: A valid merged genome fasta (UPXXXXXXXXX.fa)
    taxid: A valid genome taxonomy id
    source: The database where the genome was downloaded from
    filename: A filename to be used for the output file. If None, uses fasta_file name
    to_file: A boolean value indicating whether to write output to file or print the entries
    on screen

    returns: void
    """

    mol_type = "genomic DNA"
    previous_acc = ''

    # create an output file pointed
    output_fp = None

    if to_file is True:
        if filename is None:
            filename = os.path.basename(fasta_file).partition('.')[0]
        destination = os.path.split(fasta_file)[0]
        output_fp = open(os.path.join(destination, filename+".rfamseq"), 'w')

    args = [gf.ESL_SEQSTAT, "-a", "--dna", fasta_file]
    # open a process pipe and run esl-seqstat. Stores the result in process
    process = subprocess.Popen(args, stdout=subprocess.PIPE)

    # fetch esl-seqstat results
    result = process.communicate()

    # read in and clean-up sequence stats
    seq_stats = [x[2:] for x in result[0].strip().split('\n') if x[0] == '=']

    # process every sequence header to generate a rfamseq entry
    for sequence in seq_stats:
        seq_elements = [x.strip() for x in sequence.split(' ') if x != '']

        rfamseq_acc = seq_elements[0].split('|')[-1]
        accession = rfamseq_acc.partition('.')[0]
        version = rfamseq_acc.partition('.')[2]
        length = seq_elements[1]
        description = " ".join(seq_elements[2:])

        rfamseq_entry = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (rfamseq_acc, accession,
                                                      version, taxid,
                                                      mol_type, length,
                                                      description, previous_acc,
                                                      source)
        if to_file is True:
            output_fp.write(rfamseq_entry+'\n')
        else:
            print rfamseq_entry

    if to_file is True:
        output_fp.close()

# -----------------------------------------------------------------------------


def main(project_dir, upid_list, upid_gca_tax_file):
    """
    Main function that does some input parsing and calls
    generate_rfamseq_metadata_from_fasta to generate new entries for rfamseq table

    project_dir: The path to a project directory where the genomes are initially
    downloaded
    upid_list: A file containing a list of upids for which to generate the .rfamseq
    files
    upid_gca_tax_file: A json file generated from the latest proteomes-all.tsv
    Uniprot file

    return:
    """
    fp = open(upid_gca_tax_file, 'r')
    upid_gca_tax_dict = json.load(fp)
    fp.close()

    fp = open(upid_list, 'r')
    upids = [x.strip() for x in fp]
    fp.close()

    for upid in upids:
        subdir = os.path.join(project_dir, upid[-3:])
        updir = os.path.join(subdir, upid)
        upfasta = os.path.join(updir, upid+'.fa')

        # do some sanity checks
        if os.path.exists(upfasta):
            # defining the source of the sequences according to assembly accession
            source = "UNIPROT; ENA"
            if upid in upid_gca_tax_dict:
                if upid_gca_tax_dict[upid]["GCA"] != -1:
                    if upid_gca_tax_dict[upid]["GCA"][0:3] == "GCF":
                        source = "NIH; NCBI"
                    else:
                        source = "UNIPROT; ENA"
            else:
                print "%s not in the current Uniprot release"%upid
                continue

            generate_rfamseq_metadata_from_fasta(upfasta, upid_gca_tax_dict[upid]["TAX"],
                                    source, filename=None, to_file=True)

# -----------------------------------------------------------------------------


if __name__ == '__main__':

    project_dir = sys.argv[1]
    upid_list = sys.argv[2]
    upid_gca_tax_file = sys.argv[3]

    main(project_dir, upid_list, upid_gca_tax_file)
