#!/usr/bin/env bash

set -euo pipefail
IFS=$'\n\t'

family="$1"
rfamseq="$2"
host="$3"
port="$4"
user="$5"
password="$6"
name="$7"

sql="$(cat<<EOS
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
    AND fr.rfam_acc = "$family"
;
EOS
)"

mysql \
  --host "$host" \
  --port "$port" \
  --user "$user" \
  "-p$password" \
  --database "$name" \
  --skip-column-names \
  <<< $sql > ids

esl-sfetch -Cf $rfamseq ids
