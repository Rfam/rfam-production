import os
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
    parser.add_argument("--dest-dir", help="Destination directory to generate output to",
                        action="store", default=os.getcwd())
    parser.add_argument("--create-dump", help="Generates a JSON (.json) dump in destination directory",
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

    if args.old_rfam:
        family_accessions = extract_rfam_family_accessions(args.report)
        if args.create_dump:
            fp = open(os.path.join(args.dest_dir, "mirna_rfam_accessions_to_update_"+".json"), "w")
            try:
                json.dump(family_accessions, fp)
            except:
                print ("Unable to create json file!")
            fp.close()