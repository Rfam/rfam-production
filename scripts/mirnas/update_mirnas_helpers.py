import csv
import json

from scripts.mirnas.mirna_config import MIRNAS_CSV


def get_mirna_dict(csv_file=None):
    csv_file = MIRNAS_CSV if csv_file is None else csv_file
    with open(csv_file, mode='r') as infile:
        reader = csv.reader(infile)
        mirnas_dict = {rows[0]: {rows[1]: rows[2]} for rows in reader}
    return mirnas_dict


def get_rfam_accs(csv_file=None):
    csv_file = MIRNAS_CSV if csv_file is None else csv_file
    mirnas_dict = get_mirna_dict(csv_file)
    rfam_accs = []
    for mirna_entry in mirnas_dict.values():
        rfam_accs.append(mirna_entry.keys()[0])
    return rfam_accs


def get_mirna_ids(input_args):
    with open(input_args, 'r') as fp:
        ids_thresholds = json.load(fp)
    mirna_ids = []
    for entry in ids_thresholds:
        # entry_id = entry.keys()[0]
        # threshold = entry.values()[0]
        mirna_ids.append(entry)
    return mirna_ids
