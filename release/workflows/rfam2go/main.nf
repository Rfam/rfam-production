process BUILD {
  output:
  path('rfam2go')

  script:
  """
  perl ${params.perl_path}/export/rfam2go.pl
  """
}

process ADD_HEADER {
  input:
  path("headerless")

  output:
  path("rfam2go")

  script:
  """
  cat <<EOS headerless >rfam2go
  !version date: ${params.release.date}
  !description: A mapping of GO terms to Rfam release ${params.release.version}
  !external resource: https://rfam.org/
  !citation: Ontiveros-Palacios et al. (2025) Nucl. Acids Res. D1: D258–D267
  !contact: rfam-help@ebi.ac.uk
  !
  EOS
  """
}

process checksum {
  input:
  path(rfam2go)

  output:
  path('md5.txt')

  script:
  """
  md5sum $rfam2go > md5.txt
  """
}

workflow GENERATE_RFAM2GO {
  emit:
    rfam2go
    rfam2go_md5
  main:
    BUILD | ADD_HEADER | set { rfam2go }
    rfam2go | MD5 | set { rfam2go_md5 }
}
