process GENERATE_SEED_FILE {
  memory { acc == "RF00005" ? "20GB" : "2.5GB" }
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
  path("family*.seed")

  output:
  path('Rfam.seed'), emit: seed
  path("Rfam.seed.gz"), emit: seed_gz

  """
  find . -name 'family*.seed' | xargs -I {} cat {} > Rfam.seed
  gzip -k Rfam.seed
  """
}

workflow GENERATE_SEED {
  take:
    families
  emit:
    seed
  main:
    families \
      | GENERATE_SEED_FILE \
      | collect \
      | MERGE_SEEDS
    MERGE_SEEDS.out.seed | set { seed }
  publish:
    MERGE_SEEDS.out.seed_gz >> 'seed'
}
