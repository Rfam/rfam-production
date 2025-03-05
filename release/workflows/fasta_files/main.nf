process GENERATE_FASTA {
  tag "${acc}"
  maxForks 100

  input:
  tuple val(acc), path(rfam_seed)

  output:
  val("${acc}.fa.gz")

  """
  fasta_file_generator.py ${params.rfamseq.combined_fasta} $rfam_seed $acc ${acc}.fa
  gzip ${acc}.fa
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
    families
    rfam_seed
  emit: done
  main:
    families \
    | combine(rfam_seed) \
    | GENERATE_FASTA \
    | set { fasta }

    fasta \
    | collect \
    | COMBINE_FASTA \
    | map { 'fasta' }
    | set { done }

  publish:
    fasta >> 'fasta_files'
    COMBINED_FASTA.out >> 'fasta_files'
}
