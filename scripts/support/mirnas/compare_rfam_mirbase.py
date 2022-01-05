"""
Generate the data used in the microRNA dashboard for the following tables:
- Rfam miRNAs without matches
- Rfam non-miRNA families matching miRBase

Requires a mapping between all sequences from miRBase and the current Rfam models.
The mapping can be generated as follows:

wget https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/sequences/by-database/mirbase.fasta .
wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz .
gunzip Rfam.cm.gz
cmpress Rfam.cm
bsub -M 16000 -n 4 -o lsf.out -e lsf.err 'cmscan -o cmscan.output.txt --tblout mirbase-cmscan.tblout --cut_ga --rfam Rfam.cm mirbase.fasta'
"""

import os
import re
import subprocess
import sys

import xml.etree.ElementTree as ET
from collections import defaultdict

from utils import db_utils as db

import requests


def get_last_modified_date(rfam_acc):
    """
    Get the last modified date of a SEED file in Rfam SVN.
    """
    cmd = 'svn log --xml  https://xfamsvn.ebi.ac.uk/svn/data_repos/trunk/Families/{}/SEED'.format(rfam_acc)
    svn_output = subprocess.check_output(cmd, shell=True)
    root = ET.fromstring(svn_output)
    return root[0][1].text[:10]


def get_rfam_microrna_families_without_mirbase_hits(rfam_microrna_accs, rfam_cmscan_accs):
    """
    Compare all Rfam microRNA families from the current release with a list of
    Rfam accessions that match miRBase sequences when analysed with cmscan.
    """
    s1 = set([x['rfam_acc'] for x in rfam_microrna_accs])
    s2 = set(rfam_cmscan_accs)
    diff = s1 - s2
    print('Found {} Rfam microRNA families without miRBase hits'.format(len(diff)))
    for rfam_acc in diff:
        metadata = db.fetch_family_metadata(rfam_acc)
        line = '\t'.join([
            '=HYPERLINK("https://rfam.org/family/{0}", "{0}")'.format(rfam_acc),
            metadata['rfam_id'],
            metadata['type'],
            metadata['description'],
            get_last_modified_date(rfam_acc),
        ])
        print(line)


def get_rfam_non_microrna_families_matching_mirbase(rfam_accs, rfam_microrna_accs, rfam_cmscan_accs):
    """
    Get a list of non-microRNA families that match miRBase sequences.
    """
    s1 = set(rfam_accs)
    s2 = set([x['rfam_acc'] for x in rfam_microrna_accs])
    rfam_non_microrna_accs = s1 - s2
    diff = rfam_non_microrna_accs.intersection(rfam_cmscan_accs)
    print('Found {} Rfam non-microRNA families matching miRBase sequences'.format(len(diff)))
    return diff


def parse_cmscan_tblout(filename):
    """
    SNORA74   RF00090   URS00007E391B_9796   -          cm      103      201        1       99      +    5'    2 0.49   0.0   97.9   3.2e-26 !   Small nucleolar RNA SNORA74
    """
    data = defaultdict(list)
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = re.split(r'\s+', line)
            if re.match(r'RF\d{5}', parts[1]):
                data[parts[1]].append(parts[2])
    return data


def get_rnacentral_mirbase_xrefs():
    """
    URS000000076D	MIRBASE	MI0000019	6239	pre_miRNA
    URS00000007CB	MIRBASE	MIMAT0023447	6289	miRNA
    URS0000000EAC	MIRBASE	MIMAT0004540	10090	miRNA
    """
    filename = 'mirbase.tsv'
    data = {}
    if not os.path.exists(filename):
        cmd = 'wget https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/id_mapping/database_mappings/mirbase.tsv'
        os.system(cmd)
    with open(filename, 'r') as f:
        for line in f:
            parts = line.split('\t')
            data[parts[0]+'_'+parts[3]] = parts[2]
    return data


def get_rnacentral_description(urs_taxid):
    """
    Fetch sequence description using RNAcentral text search API.
    """
    description = ''
    url = 'https://www.ebi.ac.uk/ebisearch/ws/rest/rnacentral/entry/{}?format=json&fields=description'
    data = requests.get(url.format(urs_taxid))
    if data.status_code == 200:
        description = data.json()['entries'][0]['fields']['description'][0]
    return description


def display_non_mirna_matching_accs(non_mirna_matching_accs, cmscan_data):
    mirbase_xrefs = get_rnacentral_mirbase_xrefs()
    for rfam_acc in non_mirna_matching_accs:
        for seq in cmscan_data[rfam_acc]:
            metadata = db.fetch_family_metadata(rfam_acc)
            if seq in mirbase_xrefs:
                mirbase_url = '=HYPERLINK("https://www.mirbase.org/cgi-bin/mirna_entry.pl?acc={0}", "{0}")'.format(mirbase_xrefs[seq])
            else:
                mirbase_url = ''
            line = '\t'.join([
                '=HYPERLINK("https://rfam.org/family/{0}", "{0}")'.format(rfam_acc),
                metadata['rfam_id'],
                metadata['type'],
                metadata['description'],
                '=HYPERLINK("https://rnacentral.org/rna/{0}", "{0}")'.format(seq),
                mirbase_url,
                get_rnacentral_description(seq),
            ])
            print(line)


def main():
    if len(sys.argv) < 2 or not os.path.exists(sys.argv[1]):
        filename = os.path.join(os.path.dirname(__file__), 'mirbase-cmscan-rfam-14-6.tblout')
        print('Analysing file {}'.format(filename))
    else:
        filename = sys.argv[1]
    cmscan_data = parse_cmscan_tblout(filename)
    rfam_cmscan_accs = cmscan_data.keys()
    rfam_microrna_accs = db.fetch_mirna_families()
    rfam_accs = db.fetch_rfam_accs_sorted()
    get_rfam_microrna_families_without_mirbase_hits(rfam_microrna_accs, rfam_cmscan_accs)
    non_mirna_matching_accs = get_rfam_non_microrna_families_matching_mirbase(rfam_accs, rfam_microrna_accs, rfam_cmscan_accs)
    display_non_mirna_matching_accs(non_mirna_matching_accs, cmscan_data)


if __name__ == '__main__':
    main()
