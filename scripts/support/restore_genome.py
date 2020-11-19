import os
import sys
import json
import shutil

from scripts.export.genomes import genome_fetch as gf

# ----------------------------------------------------------------------


def redownload_genome_from_gca_report(updir, gca_report_file):
    """
    Re-downloads a genome using the GCA report file located in
    the upid directory

    return:
    """

    if not os.path.exists(updir):
        os.mkdir(updir)

    fp = open(gca_report_file, 'r')

    # parse and store accessions in a list
    accessions = [x.strip().split('\t')[0] for x in fp]

    fp.close()

    # remove GCA report header
    accessions.pop(0)

    seq_dir = os.path.join(updir, "sequences")

    # create directory or clean up old download
    if not os.path.exists(seq_dir):
        os.mkdir(seq_dir)
    else:
        shutil.rmtree(seq_dir)
        os.mkdir(seq_dir)
        os.chmod(seq_dir, 0777)

    for accession in accessions:
        gf.fetch_ena_file(accession, "fasta", seq_dir, compressed=False)

# ----------------------------------------------------------------------

def redownload_genome_from_uniprot_json(updir, upid_accession_file):
    """
    Re-downloads a genome using the upid_accessions.json file
    located in the upid directory

    return:
    """

    if not os.path.exists(updir):
        os.mkdir(updir)

    fp = open(upid_accession_file, 'r')
    acc_dict = json.load(fp)
    fp.close() 

    # parse and store accessions in a list
    accessions = acc_dict["OTHER"].values()

    seq_dir = os.path.join(updir, "sequences")

    # create directory or clean up old download
    if not os.path.exists(seq_dir):
        os.mkdir(seq_dir)
    else:
        shutil.rmtree(seq_dir)
        os.mkdir(seq_dir)
        os.chmod(seq_dir, 0777)

    for accession in accessions:
        gf.fetch_ena_file(accession, "fasta", seq_dir, compressed=False)

    if acc_dict["WGS"] != -1:
	gf.copy_wgs_set_from_ftp(acc_dict["WGS"], seq_dir)

# ----------------------------------------------------------------------

if __name__ == '__main__':

    updir = sys.argv[1]
    input_file = sys.argv[2]

    if "--ena" in sys.argv:
    	redownload_genome_from_gca_report(updir, input_file)
    elif "--uniprot" in sys.argv:
	redownload_genome_from_uniprot_json(updir, input_file)

    else:
	print "\nWrong input! Your command should look like:"
	print "\npython restore_genome.py /path/to/updir /path/to/input_file db_option"
	print "\n\n db_option: --ena|--uniprot"

