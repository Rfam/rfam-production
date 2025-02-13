process GENERATE_SEED_FILE {
  tag "${acc}"
  maxForks 50

  input:
  val(acc)

  output:
  path("${acc}.seed")

  """
  writeAnnotatedSeed.pl '${acc}'
  mv ${acc} ${acc}.seed
  """
}

process MERGE_SEEDS {
  input:
  path("*.seed")

  output:
  path('Rfam.seed')

  """
  find . -name '*.seed' | xargs -I {} cat {} >> Rfam.seed
  """
}

workflow GENERATE_SEED {
  take:
    families
  emit:
    rfam_seed
  main:
    families \
      | GENERATE_SEED_FILE \
      | collect \
      | MERGE_SEEDS \
      | set { rfam_seed }
  publish:
    rfam_seed >> 'seed'
}
