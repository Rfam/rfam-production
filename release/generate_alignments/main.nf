process BUILD_ALIGNMENT {
  tag { "$family" }

  input:
  val(family)

  output:
  path('full-alignment.sto')

  """
  rfco.pl '$family'
  cd '$family'
  rfbuild -a
  cd ..
  esl-alimerge $family/SEED $family/align > ../full-alignment.sto
  """
}

workflow GENERATE_ALIGNMENTS {
  take:
    families
  emit:
    alignments
  main:
    families | BUILD_ALIGNMENT | set { alignments }
  publish:
    alignments >> 'full_alignments'
}
