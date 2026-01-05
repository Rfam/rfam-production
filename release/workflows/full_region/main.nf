process BUILD {
  input:
  path(query)

  output:
  path("Rfam.full_region.gz")

  script:
  """
  mysql -s \
    --host=${params.db.host} \
    --port=${params.db.port} \
    --user=${params.db.user} \
    --database=${params.db.name} \
    --password=${params.db.password} \
    < $query > Rfam.full_region
  gzip Rfam.full_region
  """
}

workflow GENERATE_FULL_REGION {
  main:
    Channel.fromPath("${moduleDir}/sql/full_region.sql") | BUILD | set { full_region }
  emit:
    full_region
}
