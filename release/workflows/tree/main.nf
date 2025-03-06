process GENERATE {
  memory { acc == "RF00005" ? "20GB" : "2.5GB" }
  tag "${acc}"
  maxForks 50

  input:
  val(acc)

  output:
  val("${acc}.seed_tree")

  """
  writeAnnotatedTree.pl '${acc}'
  mv '${acc}'.taxtree '${acc}'.seed_tree
  """
}

workflow GENERATE_TREE {
  take:
    families
  emit:
    trees
  main:
    families | GENERATE | set { trees }
  publish:
    trees >> 'Rfam.seed_tree'
}
