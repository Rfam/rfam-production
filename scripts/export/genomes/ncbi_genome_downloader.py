import sys
import genome_fetch as gf

# -----------------------------------------------------------------------------


def gcf_report_parser(gcf_report):
    """
    Parses a GCF report file from NCBI and returns all the NCBI accessions

    gcf_report:

    return: a list of wgs accessions listed in the gcf report file
    """

    fp = open(gcf_report, 'r')

    accessions = [x.strip().split('\t')[4] for x in fp if x[0] != '#']

    return accessions

# -----------------------------------------------------------------------------


if __name__ == '__main__':

    gcf_report_file = sys.argv[1]
    dest_dir = sys.argv[2]

    accessions = gcf_report_parser(gcf_report_file)

    for accession in accessions:
        gf.download_fasta_from_ncbi(accession, dest_dir=dest_dir)
