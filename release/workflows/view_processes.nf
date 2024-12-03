nextflow.enable.dsl=2

process fetch_families_and_uuids {
  input:
  val(_flag)

  output:
  file('family_info')

  """
  mysql \
    --host "${params.db.host}" \
    --port "${params.db.port}" \
    --user "${params.db.user}" \
    "-p${params.db.password}" \
    --database "${params.database.name}" <<< "select rfam_acc, uuid from _post_process where status='PEND'" > family_info
  """
}

process run_rfam_view {
  tag {"$family"}
  maxForks 20
  memory { ["RF00002", "RF00005", "RF00177", "RF02542"].contains(family) ? "40GB" : "10GB" }
  errorSta ma

  input:
  tuple val(family), val(uuid)

  output:
  val('done')

  """
  ${params.perl_path}/view/rfam_family_view.pl -id $uuid -f $family family
  mysql -s `python ${params.rfamprod}/scripts/view/mysql_options.py $params.db` <<< "select rfam_acc, uuid from _post_process where status='PEND' and uuid = '$uuid'" > done
  check_empty.sh done
  """
}

workflow view_processes {
  take: start
  main:
    start \
    | fetch_families_and_uuids \
    | splitCsv(sep: "\t") \
    | run_rfam_view \
    | collect \
    | map { "Done view" } \
    | ifEmpty("No views") \
    | view \
    | set { done }
  emit: done
}

workflow {
  view_processes(Channel.of('start'))
}
