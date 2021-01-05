import json
import argparse


# ------------------------------------------------------------------------------

def parse_arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument("--report", help="Raw miRNA report in .tsv format", action="store")
    parser.add_argument("--old-rfam", help="Fetches old Rfam miRNA accessions to be updated",
                        action="store_true", default=False)
    parser.add_argument("--novel", help="Fetches all novel miRNA accessions (NOT in Rfam)",
                        action="store_true", default=False)



    return parser

# ------------------------------------------------------------------------------


def extract_rfam_family_accessions(report_file):

    fp = open(report_file, 'r')

    family_accessions = {}

    for line in fp:
        line = line.strip().split('\t')

        overlap = float(line[4])
        if overlap < 100.0:
            if line[6].upper() != "DONE":
                accession = line[3]
                if len(line.split(',')) != 2:
                    family_accessions[line[3]] = {"mirbase_id": line[0],
                                                  "threshold": float(line[1]),
                                                  "overlap": float(line[4])}

    fp.close()

    return family_accessions

# ------------------------------------------------------------------------------


if __name__ == '__main__':

    parser = parse_arguments()
    args = parser.parse_args()

    #to_commit_fp = "/hps/nobackup/production/xfam/rfam/RELEASES/14.3/miRNA_relabelled/support_code/data/new_mirnas_2020-11-17.json"
    #committed_fp = "/hps/nobackup/production/xfam/rfam/RELEASES/14.3/miRNA_relabelled/support_code/data/committed_mirnas.json"

    """
    fp = open(to_commit_fp, 'r')
    to_commit = json.load(fp)
    fp.close()

    fp = open(committed_fp, 'r')
    committed = json.load(fp)
    fp.close()

    novel = {}

    for mirna in to_commit.keys():
        if mirna not in committed:
            if mirna not in novel:
                novel[mirna] = ""
        else:
            print mirna
    print "to_commit: ", len(to_commit.keys())
    print "committed: ", len(committed.keys())
    print "novel: ", len(novel.keys())
    """

    family_accessions = extract_rfam_family_accessions(args.report)
    