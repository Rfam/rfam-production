process FETCH_UUIDS {
  input:
  val(_flag)

  output:
  file('family_info')

  """
  mysql -s \
    --host=${params.db.live.host} \
    --port=${params.db.live.port} \
    --user=${params.db.live.user} \
    --database=${params.db.live.database} \
    --password=${params.db.live.password} \
  <<< "select rfam_acc, uuid from _post_process where status='PEND'" > family_info
  """
}

process VIEW_PROCESS {
  tag { "$acc" }
  memory { ["RF00002", "RF00005", "RF00177", "RF02542"].contains(acc) ? "10GB" : "20GB" }
  maxForks 20
  when params.view_process.run

  input:
  tuple val(acc), val(uuid)

  """
  rfam_family_view.pl -id $uuid -f $acc family
  """
}

workflow RUN_VIEW_PROCESS {
  main:
    FETCH_UUIDS \
    | splitCsv(sep: "\t")
    | VIEW_PROCESS
}
