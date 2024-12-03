select rfam_acc, pdb_id, chain, pdb_start, pdb_end, bit_score, evalue_score, cm_start, cm_end, hex_colour
from pdb_full_region
where is_significant=1
order by rfam_acc;
