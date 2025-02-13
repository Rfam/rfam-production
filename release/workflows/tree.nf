process GENERATE {
  tag "${acc}"
  memory '10GB'

  input:
  val(acc)

  output:
  val("${acc}.tree")

  """
  writeAnnotatedTree.pl '${acc}'
  mv '${acc}'.taxtree '${acc}'.seed_tree
  """
}

process MERGE_TREE {
  input:
  path(trees)

  output:
  path("Rfam.seed_tree.tar.gz")

  """
  mkdir Rfam.seed_tree
  cp ${trees} Rfam.seed_tree
  tar -cf Rfam.seed_tree.tar.gz Rfam.seed_tree
  """
}

workflow GENERATE_TREE {
  take:
    families
  emit:
    tree
  main:
    families | GENERATE | collect | MERGE_TREE | set { tree }
  publish:
    tree >> 'seed_tree'
}
