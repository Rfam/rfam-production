"""
A script for checking data integrity of the Rfam-miRBase sync.
"""

from utils import db_utils as db


def verify_thresholds():
    """
    Manually export a file with manually assigned thresholds for each family,
    then compare with the thresholds of microRNA families in the database.
    Example data:
    60	RF00253
    52	RF00658
    57	RF00664
    70	RF00666
    """
    data = {}
    with open('rfam-microrna-manual-thresholds.tsv') as f:
        for line in f:
            if not 'RF' in line:
                continue
            parts = line.strip().split('\t')
            data[parts[1]] = float(parts[0])
    print('Found {} manual thresholds'.format(len(data.keys())))
    rfam_data_temp = db.fetch_mirna_families()
    rfam_data = {}
    for entry in rfam_data_temp:
        rfam_data[entry['rfam_acc']] = entry['gathering_cutoff']
    matches = 0
    for rfam_acc in data.keys():
        if rfam_acc == 'RF00273':
            continue
        if data[rfam_acc] != rfam_data[rfam_acc]:
            print('{}: {} does not match {}'.format(rfam_acc, data[rfam_acc], rfam_data[rfam_acc]))
        else:
            matches += 1
    print('Found {} matches'.format(matches))


def main():
    verify_thresholds()

if __name__ == '__main__':
    main()
