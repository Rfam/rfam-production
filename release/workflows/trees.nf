process tree_file {
  tag "${accession}"
  memory { accession == "RF00005" ? 40.GB : 10.GB }
  maxForks 10

  input:
  val(accession)

  output:
  path("${accesion}.tree")

  """
  writeAnnotatedTree.pl "${accession}"
  """
}

process tree_tarball {
  publishDir "$params.release_ftp", mode: "copy"

  input:
  path(files)

  output:
  path("Rfam.seed_tree.tar.gz")

  """
  mkdir Rfam.seed_tree
  cp ${files} Rfam.seed_tree
  tar -cf Rfam.seed_tree.tar.gz Rfam.seed_tree
  """
}

workflow trees {
  take:
    accessions
  emit:
    tree
  main:
    accessions \
    | tree_file \
    | collect \
    | tree_tarball \
    | set { tree }
}
