process BUILD {
  input:
  path(query)

  output:
  path("Rfam.pdb.gz")

  script:
  """
  mysql -s \
    --host=${params.db.host} \
    --port=${params.db.port} \
    --user=${params.db.user} \
    --database=${params.db.name} \
    --password=${params.db.password} \
    < $query > Rfam.pdb
  gzip Rfam.pdb
  """
}

workflow GENERATE_PDB {
  main:
    Channel.fromPath("${moduleDir}/sql/pdb.sql") | BUILD | set { pdb }
  emit:
    pdb
}
