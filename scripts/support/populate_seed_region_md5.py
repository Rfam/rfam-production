import sys
import hashlib
import subprocess

from config import rfam_config
from utils import seed_region as sr

ESL_PATH = rfam_config.ESL_PATH


# -------------------------------------------------------------------------

def fetch_sequence(seq_file, seq_acc, seq_start, seq_end):
    """
    Extracts a sequence from sequence file seq_file using rfamseq_acc
    and sequence start-end positions (seq_start, seq_end)
    seq_file: A sequence file in fasta format to extract a sequence from
    seq_acc: The accession of the sequence to extract
    seq_start: The starting position of the sequence/subsequence
    seq_end: The end position of the sequence/subsequence

    return: A string which corresponds to the extracted sequence
    """
    cmd = "%s %s %s/%s-%s" % (ESL_PATH,
                                 seq_file, str(seq_acc),
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

    for seq_region in seed_region_rows:
        seq_acc = seq_region[1]
        seq_start = seq_region[2]
        seq_end = seq_region[3]

        # extract sequence from sequence file
        sequence = fetch_sequence(seq_file, seq_acc, seq_start, seq_end)

        # Replace Us with Ts as done by RNAcentral
        sequence = sequence.replace('U', 'T')
        
	if sequence != '' and len(sequence) > 2:
		# convert sequence to md5
        	m = hashlib.md5()
        	m.update(sequence)
        	seq_md5 = m.hexdigest()
		# update tuple
        	rfam_acc = seq_region[0]
        	extended_tuple = (seq_md5, rfam_acc, seq_acc, seq_start, seq_end)
        	# add to new list
        	sql_data_list.append(extended_tuple)
        	extended_tuple = None

	else:
		print seq_region
    # update seed_region md5s
    #sr.update_seed_region_md5s(sql_data_list)




