process BUILD {
  input:
  path(query)

  output:
  path("Rfam.clanin")

  """
  mysql -s \
    --host=${params.db.live.host} \
    --port=${params.db.live.port} \
    --user=${params.db.live.user} \
    --database=${params.db.live.database} \
    --password=${params.db.live.password} \
    < "${query}" > Rfam.clanin
  """
}

workflow GENERATE_CLANIN {
  emit:
    clanin
  main:
    channel.fromPath('sql/clanin.sql') \
    | GENERATE \
    | set { clanin }
}
