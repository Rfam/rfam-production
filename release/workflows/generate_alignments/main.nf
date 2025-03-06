process BUILD_ALIGNMENT {
  tag { "$acc" }

  input:
  val(acc)

  output:
  path('full-alignment.sto')

  """
  rfco.pl '$acc'
  cd '$acc'
  rfbuild -a
  cd ..
  esl-alimerge $acc/SEED $acc/align > ../full-alignment.sto
  """
}

workflow GENERATE_ALIGNMENTS {
  take:
    fasta
    seed
  emit:
    alignments
  main:
    families | BUILD_ALIGNMENT | set { alignments }
  publish:
    alignments >> 'full_alignments'
}
