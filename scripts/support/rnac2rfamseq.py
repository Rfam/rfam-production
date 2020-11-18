import os
import sys
import subprocess
import requests
import argparse
import csv

from config import rfam_config as rc

"""
Metadata queries:

ZWD
---

select
rna.len,
xref.upi,
xref.taxid,
acc.accession,
acc.description,
acc.feature_name,
acc.ncrna_class
from xref
join rnc_accessions acc on acc.accession = xref.ac
join rna on rna.upi = xref.upi
where
  acc.database = 'ZWD'
  and xref.deleted = 'N'


PDBe
----

select
r.len,
r.upi,
pre.taxid,
'accession' as accession,
pre.description,
pre.rna_type as feature_name,
pre.rna_type
 from rna r
 join rnc_rna_precomputed pre on r.upi = pre.upi
 where
(pre.databases like '%PDBe%' or  pre.databases like 'PDBe%' or  pre.databases like '%PDBe')
 and pre.is_active  = true
 and taxid is not null ;
"""

INSDC_MOL_TYPES = ['protein', 'genomic DNA', 'DNA', 'ss-DNA', 'RNA', 'genomic RNA',
                   'ds-RNA', 'ss-cRNA', 'ss-RNA', 'mRNA', 'tRNA', 'rRNA', 'snoRNA',
                   'snRNA', 'scRNA', 'pre-RNA', 'other RNA', 'other DNA',
                   'unassigned DNA', 'unassigned RNA', 'viral cRNA', 'cRNA',
                   'transcribed RNA', 'ncRNA', 'ribozyme', 'antisense_RNA', 'other']

# ----------------------------------------------------------------------------------


def load_rnacentral_metadata_to_dict(rnac_metadata_file, to_file=False, destdir=None):
    """
    Loads RNAcentral metadata to a dictionary. The expected format is:
    URS_acc\tsource\tpre_acc/start-end\tncbi_id\tmol_type

    rnac_metadata_file: A tab delimited file with Rfam compatible metadata
    extracted from RNAcentral
    to_file: If True, it dumps the output dictionary in a json file
    destdir: The path to a destination directory where any output will be generated

    return: A dictionary with all metadata values loaded from rnacentral tsv file
    """

    # open a new file handle
    fp = open(rnac_metadata_file, 'r')

    rnac_metadata_dict = {}
    # loop over all lines in the rnacentral tsv file
    for rnac_line in fp:
        rnac_line_contents = rnac_line.strip().split('\t')
        tax_id = rnac_line_contents[3]

        # append taxid to URS to create RNAcentral accession
        rnac_acc = rnac_line_contents[0] + '_' + tax_id
        # all URS entries must be unique, count if any skipped
        if rnac_acc not in rnac_metadata_dict:
            rnac_metadata_dict[rnac_acc] = {}
            rnac_metadata_dict[rnac_acc]["source"] = rnac_line_contents[1]
            temp = rnac_line_contents[2].partition('/')

            previous_acc = temp[0]
            """ Revise this section
            if temp[2] != '':
                region_coords = temp[2].split('-')
                seq_start = region_coords[0]
                seq_end = region_coords[1]
            else:
                seq_start = 1
                seq_end = seq_length
            """
            rnac_metadata_dict[rnac_acc]["previous_acc"] = previous_acc
            rnac_metadata_dict[rnac_acc]["tax_id"] = tax_id
            rnac_metadata_dict[rnac_acc]["mol_type"] = rnac_line_contents[4]

    # close file handle
    fp.close()

    # create a json dump
    if to_file is True:
        import json
        filename = os.path.basename(rnac_metadata_file).partition('.')[0] + '.json'

        if destdir is None:
            destdir = os.path.split(rnac_metadata_file)[0]

        json_dump_file = os.path.join(os.path.join(destdir, filename))

        fp = open(json_dump_file, 'w')
        json.dump(rnac_metadata_dict, fp)

    return rnac_metadata_dict

# ----------------------------------------------------------------------------------


def generate_sequence_stats(fasta_file):
    """
    Uses esl-sfetch to generate stats for every sequence in the provided fasta file
    and loads that information in a dictionary

    fasta_file: A sequence file in fasta format

    return: A dictionary with sequence statistics
    """

    sequence_stats = {}
    # create argument list required by subprocess PIPE
    args = [rc.ESL_SEQSTAT, "-a", "--dna", fasta_file]

    # open a process pipe and run esl-seqstat. Stores the result in process
    process = subprocess.Popen(args, stdout=subprocess.PIPE)

    # fetch esl-seqstat results
    result = process.communicate()

    # read in and clean-up sequence stats
    seq_stats = [x[2:] for x in result[0].strip().split('\n') if x[0] == '=']

    # process every sequence header to generate a rfamseq entry
    for sequence in seq_stats:
        seq_elements = [x.strip() for x in sequence.split(' ') if x != '']

        rfamseq_acc = seq_elements[0]

        if rfamseq_acc not in sequence_stats:
            sequence_stats[rfamseq_acc] = {}
            sequence_stats[rfamseq_acc]["length"] = seq_elements[1]
            sequence_stats[rfamseq_acc]["desc"] = ' '.join(seq_elements[2:])

    return sequence_stats

# ----------------------------------------------------------------------------------


def generate_rfamseq_tsv_dump(fasta_file, rnac_metadata_file, to_file = False):
    """
    Parses the sequence file and the rnacentral metadata file and generates an
    rfamseq dump for the sequences to be imported to rfam_live

    fasta_file: The path to the input fasta file as generated from RNAcentral
    rnac_metadata_file: The path to the RNAcentral metadata file in tabular
    format

    return: void
    """

    rnac_metadata = load_rnacentral_metadata_to_dict(rnac_metadata_file)
    sequence_stats = generate_sequence_stats(fasta_file)

    for urs_id in rnac_metadata.keys():
        rfamseq_acc = urs_id
        accession = urs_id
        version = '0'
        taxid = rnac_metadata[urs_id]["tax_id"]
        mol_type = rnac_metadata[urs_id]["mol_type"]
        length = sequence_stats[urs_id]["length"]
        description = sequence_stats[urs_id]["desc"]
        previous_acc = rnac_metadata[urs_id]["previous_acc"]
        source = rnac_metadata[urs_id]["source"]


        rfamseq_entry = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (rfamseq_acc, accession,
                                                          version, taxid,
                                                          mol_type, length,
                                                          description, previous_acc,
                                                          source)

        # TO DO provide the option to generate a file instead?

        print rfamseq_entry


# ----------------------------------------------------------------------------------


def fetch_rfamseq_metadata_from_rnacentral(urs_accession):
    """

    urs_accession:
    return:
    """

    rnacentral_url = 'https://rnacentral.org/api/v1/rna'
    response = requests.get(rnacentral_url, params={'urs': urs_accession})

    data = response.json()

    if data['count'] > 0:
        return data['results'][0]


# ----------------------------------------------------------------------------------

def rnacentral_csv_2_rfamseq(rnac_csv, source, delimiter='\t'):
    """

    rnac_csv:

    return:
    """

    # 53,URS0000BFE20F,12908,ZWD:JCVI_SCAF_1096626860522/137-85,unclassified sequences glnA RNA,ncRNA,ncRNA
    csv_in = open(rnac_csv, 'r')

    for line in csv_in:
        line = line.strip().split(delimiter)
        rfamseq_acc = line[1] + '_' + line[2]
        version = "000000"
        ncbi_id = line[2]

        mol_type = ''
        if line[6] not in INSDC_MOL_TYPES:
            mol_type = "other"
        else:
            mol_type = line[6]

        length = line[0]
        description = line[4].replace('\"', '')

        if line[3] == "accession":
            prev_acc = ''
        else:
            prev_acc = line[3]

        new_entry = '\t'.join([rfamseq_acc, rfamseq_acc, version, ncbi_id, mol_type,
                           length, description, prev_acc, source])

        print new_entry

    csv_in.close()

# --------------------------------------------------------------------------------


def rnacentral_2_rfamseq(fasta_file):
    """

    fasta_file:
    return:
    """

    fasta_fp = open(fasta_file, 'r')
    rfamseq_entry = []

    for line in fasta_fp:
        if line[0] == '>':
            urs_accession = line[1:].strip()

            urs_metadata = fetch_rfamseq_metadata_from_rnacentral(urs_accession)

            rfamseq_entry.append(urs_accession)
            rfamseq_entry.append(urs_accession)
            rfamseq_entry.append("000000")
            rfamseq_entry.append(urs_accession.partition('_')[2])
            rfamseq_entry.append(urs_metadata["rna_type"])
            rfamseq_entry.append(str(urs_metadata["length"]))
            rfamseq_entry.append(urs_metadata["description"])
            rfamseq_entry.append("ZWD")

            print ('\t'.join(rfamseq_entry))

            rfamseq_entry = []

    fasta_fp.close()

# ----------------------------------------------------------------------------------


def parse_arguments():
    """
    Basic argument parsing using python's argparse

    return: void
    """

    parser = argparse.ArgumentParser(description="Script to convert RNAcentral sequence metadata to Rfamseq")

    parser.add_argument("--rnac-csv",
                        help="A valid csv file with RNAcentral sequence metadata", action="store")
    parser.add_argument("--source",
                        help="The source/database the metadata was obtained from (e.g. ZWD)", action="store")

    return parser

# ----------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = parse_arguments()
    args = parser.parse_args()

    rnacentral_csv_2_rfamseq(args.rnac_csv, args.source)

    """

    fasta_file = sys.argv[1]
    rnac_metadata_file = sys.argv[2]

    generate_rfamseq_tsv_dump(fasta_file, rnac_metadata_file, to_file=False)

    """