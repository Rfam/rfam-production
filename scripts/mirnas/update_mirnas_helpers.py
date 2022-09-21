import csv


def get_mirna_dict(input_args):
    with open(input_args, mode='r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        mirnas_dict = {rows[0]: {rows[1]: rows[2]} for rows in reader}
    return mirnas_dict


def get_rfam_accs(input_args):
    mirnas_dict = get_mirna_dict(input_args)
    rfam_accs = []
    for mirna_entry in mirnas_dict.values():
        rfam_accs.append(mirna_entry.keys()[0])
    return rfam_accs


def get_mirna_ids(input_args):
    mirna_ids = []
    with open(input_args, mode='r') as infile:
        tsv_file = csv.reader(infile, delimiter='\t')
        for line in tsv_file:
            mirna_id = line.split('\t')
            mirna_ids.append(mirna_id)
    return mirna_ids
