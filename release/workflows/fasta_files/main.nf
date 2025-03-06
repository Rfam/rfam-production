process GENERATE_FASTA {
  tag "${acc}"
  maxForks 50

  input:
  tuple val(acc), path(rfam_seed)

  output:
  val("${acc}.fa.gz")

  """
  {
    esl-reformat fasta '${family_seed}'
    export-full.sh \
      '${acc}' \
      '${params.rfamseq.directory}/${params.rfamseq.combined_fasta}' \
      '${params.db.host}' \
      '${params.db.port}' \
      '${params.db.user}' \
      '${params.db.password}' \
      '${params.db.database}'
  } > '${acc}.unsorted.fa'
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
  main:
    seed_alignments | GENERATE_FASTA | set { fasta }

    fasta | collect | COMBINE_FASTA | set { all_fasta }
  publish:
    fasta >> 'fasta_files'
    all_fasta >> 'fasta_files'
}
