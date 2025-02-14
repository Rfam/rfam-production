process FETCH_FAMILIES {
  output:
  path('families')

  """
  mysql -s \
    --host ${params.db.live.host} \
    --port ${params.db.live.port} \
    --user ${params.db.live.user} \
    --password ${params.db.live.password} \
    <<< "select rfam_acc from family" > families
  """
}
