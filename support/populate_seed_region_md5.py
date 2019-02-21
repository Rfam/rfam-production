import hashlib
import os
import subprocess
import sys

from utils import seed_region as sr

# -------------------------------------------------------------------------

def fetch_sequence_from_seed(rfam_seed_file, seq_acc, seq_start, seq_end):
    """
    Extracts a sequence from sequence file seq_file using rfamseq_acc
    and sequence start-end positions (seq_start, seq_end)
    rfam_seed_file: A sequence file in fasta format to extract a sequence from
    seq_acc: The accession of the sequence to extract
    seq_start: The starting position of the sequence/subsequence
    seq_end: The end position of the sequence/subsequence

    return: A string which corresponds to the extracted sequence
    """
    cmd = "esl-sfetch %s %s/%s-%s" % (seq_file, str(seq_acc),
                              str(seq_start), str(seq_end))

    proc = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE)

    seq = proc.communicate()[0]

    # get sequence
    sequence = ''
    seq_bits = seq.split('\n')[1:]
    sequence = sequence.join(seq_bits)

    return sequence


# -------------------------------------------------------------------------

if __name__ == '__main__':

    seq_file = sys.argv[1]

    seed_region_rows = sr.fetch_seed_regions()

    sql_data_list = []

    fp_err = None

    for seq_region in seed_region_rows:
        seq_acc = seq_region[1]
        seq_start = seq_region[2]
        seq_end = seq_region[3]
        rfam_acc = seq_region[0]

        # extract sequence from sequence file
        sequence = fetch_sequence_from_seed(seq_file, seq_acc, seq_start, seq_end)

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
                fp_err = open(os.path.join(os.path.split(seq_file)[0], "seed_md5.err"), 'w')

            fp_err.write("Error extracting sequence %s at %s-%s of %s\n" % (seq_acc,
                                                                            seq_start,
                                                                            seq_end,
                                                                            rfam_acc))

    # update/populate table
    sr.update_seed_region_md5s(sql_data_list)
    # close error file if opened
    if fp_err is not None:
        fp_err.close()
