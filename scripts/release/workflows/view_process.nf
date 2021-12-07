nextflow.enable.dsl=2

process fetch_families_and_uuids {
  output:
  file('family_info')

  """
  mysql `python ${params.rfamprod}/scripts/view/mysql_options.py $params.db` <<< "select rfam_acc, uuid from _post_process where status='PEND'" > family_info
  """
}

process run_rfam_view {
  tag {"$family"}
  memory {
    ["RF00002", "RF00005", "RF00177", "RF02542"].contains(family) ? 10000 : 20000
  }

  clusterOptions '-g /rfam_view'

  input:
  tuple val(family), val(uuid)

  """
  ${params.perl_path}/view/rfam_family_view.pl -id $uuid -f $family family
  mysql `python ${params.rfamprod}/scripts/view/mysql_options.py $params.db` <<< "select rfam_acc, uuid from _post_process where status='PEND' and uuid = '$uuid'" > done
  bash $baseDir/check_empty.sh done
  """
}

workflow view_process {
  main:
    fetch_families_and_uuids \
    | splitCsv(sep: "\t") \
    | run_rfam_view
}

workflow {
  view_process()
}
