process GENERATE_FASTA {
  tag "${acc}"
  maxForks 50

  input:
  tuple val(acc), path(rfam_seed)

  output:
  val("${acc}.fa.gz")

  script:
  sql = """SELECT
    CONCAT(fr.rfamseq_acc, '/', fr.seq_start, '-', fr.seq_end),
    fr.seq_start,
    fr.seq_end,
    fr.rfamseq_acc,
    rf.description
FROM full_region fr
JOIN rfamseq rf USING (rfamseq_acc)
WHERE
    fr.is_significant = 1
    AND fr.rfam_acc = '${acc}';"""

  """
  mysql \
    --host='${params.db.live.host}' \
    --port='${params.db.live.port}' \
    --user='${params.db.live.user}' \
    -p'${params.db.live.password}' \
    --database='${params.db.live.database}' \
    --skip-column-names <<< "$sql" > ids

  esl-reformat fasta '${rfam_seed}'                                                  > '${acc}.unsorted.fa'
  esl-sfetch -Cf '${params.rfamseq.directory}/${params.rfamseq.combined_fasta}' ids >> '${acc}.unsorted.fa'
  seqkit rmdup '${acc}.unsorted.fa' > '${acc}.fa'
  gzip '${acc}.fa'
  """
}

process COMBINE_FASTA {
  input:
  path("family*.fa.gz")

  output:
  path('Rfam.fa.gz')

  """
  find . -name 'family*.fa.gz' | xargs -I {} zcat {} > Rfam.fa.unsorted
  seqkit rmdup Rfam.fa.unsorted > Rfam.fa
  gzip Rfam.fa
  """
}

workflow GENERATE_FASTA_FILES {
  take:
    seed_alignments
  emit:
    fasta
    all_fasta
  main:
    seed_alignments | GENERATE_FASTA | set { fasta }
    fasta | collect | COMBINE_FASTA | set { all_fasta }
}
