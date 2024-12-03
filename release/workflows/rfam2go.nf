nextflow.enable.dsl = 2

process export {
  publishDir "$params.release_ftp/rfam2go", mode: "copy"

  input:
  path(_flag)

  output:
  tuple path('rfam2go'), path('md5.txt')

  """
  perl $params.perl_path/export/rfam2go.pl
  add_header_rfam2go.py --input rfam2go --rfam-version $params.releasex > with-header
  mv with-header rfam2go
  md5sum rfam2go > md5.txt
  """
}

workflow rfam2go {
  take: start
  emit: done
  main:
    start | export | set { done }
}

workflow {
  rfam2go(Channel.of('start'))
}
