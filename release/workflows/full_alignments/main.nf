process BUILD_ALIGNMENT {
  tag { "$acc" }

  input:
  tuple val(acc), path(seed_alignment)

  output:
  path('full-alignment.sto')

  """
  rfco.pl '$acc'
  cd '$acc'
  rfbuild -a
  cd ..
  esl-alimerge '${seed_alignment}' '$acc/align' > full-alignment.sto
  """
}

workflow GENERATE_FULL_ALIGNMENTS {
  take:
    seed_alignments
  emit:
    full_alignments
  main:
    families | BUILD_ALIGNMENT | set { full_alignments }
}
