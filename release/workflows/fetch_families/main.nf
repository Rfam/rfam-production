process FETCH_FAMILIES {
  output:
  path('families')

  """
  mysql -s \
    --host=${params.db.live.host} \
    --port=${params.db.live.port} \
    --user=${params.db.live.user} \
    --database=${params.db.live.database} \
    --password=${params.db.live.password} \
    <<< "select rfam_acc from family" > families
  """
}
