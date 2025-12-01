process BUILD_ALIGNMENT {
  tag { "$acc" }
  maxForks 50

  input:
  val(acc)

  output:
  path("${acc}.sto")

  script:
  """
  rfco.pl '$acc'
  cd '$acc'
  rfbuild -a
  cd ..
  mv $acc/align "${acc}.sto"
  """
}

workflow GENERATE_FULL_ALIGNMENTS {
  take:
  accessions
  emit:
  full_alignments
  main:
  accessions | BUILD_ALIGNMENT | set { full_alignments }
}
