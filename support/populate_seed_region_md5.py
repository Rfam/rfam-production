import hashlib
import os
import subprocess
import sys

from utils import seed_region as sr
from utils import db_utils as db

# -------------------------------------------------------------------------

def fetch_sequence(seq_file, seq_acc, seq_start, seq_end, type='seed'):
    """
    Extracts a sequence from sequence file seq_file using rfamseq_acc
    and sequence start-end positions (seq_start, seq_end)
    rfam_seed_file: A sequence file in fasta format to extract a sequence from
    seq_acc: The accession of the sequence to extract
    seq_start: The starting position of the sequence/subsequence
    seq_end: The end position of the sequence/subsequence

    return: A string which corresponds to the extracted sequence
    """
    cmd = ''
    if type == 'seed':
        cmd = "esl-sfetch %s %s/%s-%s" % (seq_file, str(seq_acc),
                              str(seq_start), str(seq_end))

    elif type == 'full':
        cmd = "esl-sfetch -c %s..%s %s -%s" % (str(seq_start), str(seq_end),
                                               seq_file, str(seq_acc))
    proc = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE)

    seq = proc.communicate()[0]

    # get sequence
    sequence = ''
    seq_bits = seq.split('\n')[1:]
    sequence = sequence.join(seq_bits)

    return sequence

# -------------------------------------------------------------------------

def generate_md5s_and_populate_table(seq_file, type = "seed", dest_dir = None):
    """
    Fetches seed of full regions from the database and populates/updates the
    appropriate table with the md5s of the corresponding ncRNA sequences.

    seq_file: This can be the Rfam.seed file (seed) or an Rfamseq file (full)
    type: one of seed/full

    return: void
    """

    if dest_dir is None:
        dest_dir = os.path.split(seq_file)[0]


    region_rows = None

    if type == "seed":
        region_rows = sr.fetch_seed_regions()
    elif type == "full":
        region_rows = db.fetch_metagenomic_regions()

    sql_data_list = []

    fp_err = None

    for seq_region in region_rows:
        seq_acc = seq_region[1]
        seq_start = seq_region[2]
        seq_end = seq_region[3]
        rfam_acc = seq_region[0]

        # extract sequence from sequence file
        sequence = fetch_sequence(seq_file, seq_acc, seq_start, seq_end, type)

        # Replace Ts with Us as done by RNAcentral
        # Note: there shouldn't be any in the seed alignments
        sequence = sequence.replace('T', 'U')

        if sequence != '' and len(sequence) > 2:
            # convert sequence to md5
            m = hashlib.md5()
            m.update(sequence)
            seq_md5 = m.hexdigest()
            # update tuple
            extended_tuple = (seq_md5, rfam_acc, seq_acc, seq_start, seq_end)
            # add to new list
            sql_data_list.append(extended_tuple)
            extended_tuple = None

        else:
            # check if first error
            if fp_err is None:
                fp_err = open(os.path.join(dest_dir, "md5.err"), 'w')

            fp_err.write("Error extracting sequence %s at %s-%s of %s\n" % (seq_acc,
                                                                            seq_start,
                                                                            seq_end,
                                                                            rfam_acc))

    # update/populate table
    if type == "seed":
        sr.update_seed_region_md5s(sql_data_list)
    else:
        db.update_metagenomic_region_md5s(sql_data_list)

    # close error file if opened
    if fp_err is not None:
        fp_err.close()


# -------------------------------------------------------------------------
if __name__ == '__main__':

    seq_file = sys.argv[1]
    type_flag = sys.argv[2]

    if type_flag.lower() == "--seed":
        type = "seed"
    elif type_flag.lower() == "--full":
        type = "full"
    else:
        print "\nWrong selection of type: %s\n" % type_flag

    generate_md5s_and_populate_table(seq_file, type="seed", dest_dir=None)



