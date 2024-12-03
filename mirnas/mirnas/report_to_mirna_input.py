import argparse
import json
import os
from datetime import date

# -------------------------------------------------------------------------------


def extract_new_mirnas_from_report(report_tsv, type='new'):
    """
    """

    new_mirnas = {}

    fp = open(report_tsv, 'r')
    count = 0
    for line in fp:
        line = line.strip().split('\t')
        # check if candidate mirna is a new family
        if line[6].lower() == "new family":
            # skip families requiring review
            if line[1] != '?' and line[1] != '' and line[2] != "1_SEED":
                if line[0] not in new_mirnas:
                    print line
                    new_mirnas[line[0]] = line[1]
        elif line[6].lower() == 'done':
            count += 1
    fp.close()

    return new_mirnas


# -------------------------------------------------------------------------------

def extract_rfam_family_accessions(report_file):

    fp = open(report_file, 'r')

    accession_map = {}

    for line in fp:
        line = line.strip().split('\t')

        overlap = float(line[4])
        if overlap <= 100.0:
            # converts to upper to ensure labels match the constant
            if line[6].strip().upper() == "UPDATE SEED":

                rfam_acc = line[3].strip()

                rfam_acc_list = []
                if rfam_acc.find(',') == -1:
                    rfam_acc_list = rfam_acc.split(',')
                else:
                    rfam_acc_list = [rfam_acc]

                threshold = 0.0
                if line[1] != '':
                    threshold = float(line[1])

                # trim any whitespace characters
                mirbase_id = line[0].strip()


                accession_map[mirbase_id] = {"rfam_acc": rfam_acc_list,
                                              "threshold": threshold,
                                              "overlap": float(line[4])}

    fp.close()

    return accession_map


# -------------------------------------------------------------------------------

def parse_arguments():
    """
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("--report", help="miRNA report in .tsv format", action="store")
    parser.add_argument("--dest-dir", help="Desctination directory", action="store", default=os.getcwd())
    parser.add_argument("--old-rfam", help="Fetches old Rfam miRNA accessions to be updated",
                        action="store_true", default=False)
    parser.add_argument("--create-dump", help="Generates a JSON (.json) dump in destination directory",
                        action="store_true", default=False)

    return parser


# -------------------------------------------------------------------------------

if __name__ == '__main__':

    parser = parse_arguments()
    args = parser.parse_args()

    accessions = None
    if not args.old_rfam:
        new_mirnas = extract_new_mirnas_from_report(args.report, type='new')
        accessions = new_mirnas

    else:
        accessions = extract_rfam_family_accessions(args.report)

    if args.create_dump:
        filename = "new_mirnas_"

        if args.old_rfam:
            filename = "mirna_families_to_update_"
        fp_out = open(os.path.join(args.dest_dir, filename + str(date.today()) + ".json"), 'w')
        json.dump(accessions, fp_out)
        fp_out.close()
