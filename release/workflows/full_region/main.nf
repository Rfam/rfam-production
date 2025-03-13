process BUILD {
  input:
  path(query)

  output:
  path("Rfam.full_region.gz")

  """
  mysql -s \
    --host=${params.db.live.host} \
    --port=${params.db.live.port} \
    --user=${params.db.live.user} \
    --database=${params.db.live.database} \
    --password=${params.db.live.password} \
    < $query > Rfam.full_region
  gzip Rfam.full_region
  """
}

workflow GENERATE_FULL_REGION {
  emit:
    full_region
  main:
    channel.fromPath("${moduleDir}/sql/full_region.sql") | BUILD | set { full_region }
}
