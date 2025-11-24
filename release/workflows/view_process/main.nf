process FETCH_UUIDS {
  input:
  val(_flag)

  output:
  file('family_info')

  """
  mysql -s \
    --host=${params.db.host} \
    --port=${params.db.port} \
    --user=${params.db.user} \
    --database=${params.db.name} \
    --password=${params.db.password} \
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
    FETCH_UUIDS(true) \
    | splitCsv(sep: "\t")
    | VIEW_PROCESS
}
