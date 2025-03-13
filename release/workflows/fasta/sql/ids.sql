SELECT
    CONCAT(fr.rfamseq_acc, '/', fr.seq_start, '-', fr.seq_end),
    fr.seq_start,
    fr.seq_end,
    fr.rfamseq_acc,
    rf.description
FROM full_region fr
JOIN rfamseq rf USING (rfamseq_acc)
WHERE
    fr.is_significant = 1
    AND fr.rfam_acc = '${acc}'
;
