import os
import hashlib
import requests

# ---------------------------------------------------------------

def sequence_to_md5(sequence):
    """
    Converts a sequence to an md5 hash after replacing Us with
    Ts

    sequence: A valid RNA/DNA sequence

    return: MD5 hash of the sequence
    """

    md5_converter = hashlib.md5()
    # convert to DNA
    sequence = sequence.replace('U', 'T')
    md5_converter.update(sequence.encode('utf-8'))
    sequence_md5 = md5_converter.hexdigest()

    return sequence_md5

# ---------------------------------------------------------------


def generate_seed_id_from_RNAcentral(sequence):
    """
    Generates a seed accession based on a sequence mad5 match in RNAcentral

    sequence: A valid DNA/RNA sequence

    return: Returns RNAcentral id, otherwise returns None
    """

    sequence_md5 = sequence_to_md5(sequence)

    rnacentral_url = 'https://rnacentral.org/api/v1/rna'
    response = requests.get(rnacentral_url, params={'md5': sequence_md5})

    data = response.json()

    if data['count'] > 0:
        return data['results'][0]['rnacentral_id'] + "/1-" + str(data['results'][0]['length'])

    return None

# ---------------------------------------------------------------


def seed_to_dict(seed):
    """
    Relabels the accessions of a SEED alignment using RNAcentral
    identifiers. This is done by matching the seed sequences, with
    sequences existing in RNAcentral using md5 hashing.

    seed: A reformatted seed in Pfam format
    dest_dir: The path to the destination directory. None by default

    return: The path to the relabelled SEED alignement
    """
    rnac_miRNA_mappings = {}
    miRNA_mappings = {}

    filename = os.path.split(seed)[1].partition('.')[0]

    seed_fp = open(seed, 'r')

    for line in seed_fp:
        # check if this is an actual sequence line
        if line[0] != '#' and len(line) > 1 and line[0:2] != '//':
            line_elements = [x for x in line.strip().split(' ') if x != '']
            miRNA_mappings[line_elements[0]] = line_elements[1]

            sequence = line_elements[1].replace('.', '').replace('t', 'u').replace('T', 'U').upper()
            #md5_seq = sequence_to_md5(sequence)

            rnac_accession = generate_seed_id_from_RNAcentral(sequence)

            if rnac_accession is not None:
                rnac_miRNA_mappings[rnac_accession] = sequence

    seed_fp.close()

    return miRNA_mappings, rnac_miRNA_mappings

# ---------------------------------------------------------------


if __name__ == '__main__':

    seed = "/Users/ikalvari/Desktop/releases/14.3/all_seeds.stk"
    miRNA_mappings, rnac_miRNA_mappings = seed_to_dict(seed)



    print len(miRNA_mappings.keys())
    print len(rnac_miRNA_mappings.keys())

