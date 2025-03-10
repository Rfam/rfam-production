process BUILD {
  input:
  path(query)

  output:
  path("Rfam.pdb.gz")

  """
  mysql -s \
    --host=${params.db.live.host} \
    --port=${params.db.live.port} \
    --user=${params.db.live.user} \
    --database=${params.db.live.database} \
    --password=${params.db.live.password} \
    < $query > Rfam.pdb
  gzip Rfam.pdb
  """
}

workflow GENERATE_PDB {
  emit:
    pdb
  main:
    channel.fromPath('sql/pdb.sql') \
    | BUILD \
    | set { pdb }
}
