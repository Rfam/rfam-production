select rfam_acc, rfamseq_acc, seq_start, seq_end, bit_score, evalue_score, cm_start, cm_end, truncated, type
from full_region
where is_significant=1
order by rfam_acc;
