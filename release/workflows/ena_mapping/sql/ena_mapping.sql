select
    accession,
    rfam_acc,
    rfam_id
from full_region
join family using (rfam_acc)
join rfamseq using (rfamseq_acc)
where
    full_region.is_significant
