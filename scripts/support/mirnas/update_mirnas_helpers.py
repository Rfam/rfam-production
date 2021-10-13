import csv

MIRNAS_CSV = "rfam-production/mirnas_sample.csv"
UPDATE_DIR = "/hps/nobackup/production/xfam/rfam/RELEASES/14.3/miRNA_relabelled/update_old_rfam_mirnas"
MEMORY = 8000
CPU = 4
LSF_GROUP = "/rfam_srch"


def get_data_from_csv(csv_file=None):
    csv_file = MIRNAS_CSV if csv_file is None else csv_file
    with open(csv_file, mode='r') as infile:
        reader = csv.reader(infile)
        mirnas_dict = {rows[0]: {rows[1]: rows[2]} for rows in reader}
    return mirnas_dict


def get_rfam_accs(csv_file=None):
    csv_file = MIRNAS_CSV if csv_file is None else csv_file
    mirnas_dict = get_data_from_csv(csv_file)
    rfam_accs = []
    for mirna_entry in mirnas_dict.values():
        rfam_accs.append(mirna_entry.keys()[0])
    return rfam_accs
