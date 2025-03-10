process GENERATE_SEED_FILE {
  memory { acc == "RF00005" ? "20GB" : "2.5GB" }
  tag "${acc}"
  maxForks 50

  input:
  val(acc)

  output:
  tuple val(acc), path("${acc}.seed")

  """
  writeAnnotatedSeed.pl '${acc}'
  mv ${acc} ${acc}.seed
  """
}

process MERGE_SEEDS {
  input:
  path("family*.seed")

  output:
  path("Rfam.seed.gz"), emit: seed_gz

  """
  find . -name 'family*.seed' | xargs -I {} cat {} > Rfam.seed
  gzip Rfam.seed
  """
}

workflow GENERATE_SEED {
  take:
    families
  emit:
    seeds
    seed_gz
  main:
    families | GENERATE_SEED_FILE | set { seeds }
    seeds | map { it[1] } | collect | MERGE_SEEDS
    MERGE_SEEDS.out.seed_gz | set { seed_gz }
}
