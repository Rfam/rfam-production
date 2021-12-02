nextflow.enable.dsl=2

process fetch_families_and_uuids {
  output:
  file('family_info')

  """
  mysql $params.mysql_options <<< "select rfam_acc, uuid from _post_process where status='PEND'" > family_info
  """
}

process run_rfam_view {
  label "$family"
  memory {
    ["RF00002", "RF00005", "RF00177", "RF02542"].contains(family) ? 10000 : 20000
  }

  queueOptions '-g /rfam_view'

  input:
  tuple val(family), val(uuid)

  """
  $params.rfam_family_view -id $uuid -f $family family
  mysql $params.mysql_options <<< "select rfam_acc, uuid from _post_process where status='PEND' and uuid = '$uuid'" > done
  test -n done
  """
}

workflow view_process {
  main:
    fetch_families_and_uuids \
    | splitCsv(sep: "\t") \
    | run_rfam_view
}
