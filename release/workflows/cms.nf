#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process family_cm {
  tag "${accession}"
  memory { acc == "RF00005" ? 20.GB : 5.GB }
  maxForks 20

  input:
  val(accession)

  output:
  path("${accession}.cm")

  """
  writeAnnotatedCM.pl $accession
  mv ${accession}.CM initial.cm
  grep -v DESC initial.cm > ${accession}.nodesc.cm
  perl $params.perl_path/jiffies/seed-desc-to-cm.pl $params.release_ftp/seed/Rfam.seed Rfam.nodesc.cm > ${accession}.cm
  """
}

process tar_cm_files {
  stageInMode 'copy'
  publishDir "${params.release_ftp}", mode: "copy"

  input:
  path(cm)

  output:
  path('Rfam.tar.gz')

  // Note this assumes all input $cm are named RF*.cm
  """
  tar -czvf Rfam.tar.gz RF*.cm
  """
}

process merge_cm {

  input:
  path(cm)

  output:
  path('Rfam.cm.gz')

  """
  cat $cm > Rfam.cm
  """
}

process validate_cms {
  publishDir "${params.release_ftp}", mode: "copy"

  input:
  tuple path(families), path('Rfam.cm')

  output:
  path('Rfam.cm.gz')

  """
  cmstat Rfam.cm > cmstat_file.txt
  cm_check.py --accessions $families --cm-file Rfam.cm --stat-file cmstat_file.txt
  gzip Rfam.cm
  """
}

workflow cms {
  take:
    accession_file
    families
  main:
    families \
    | family_cm \
    | collect \
    | set { cm_list }

    cm_list | tar_cm_files | set { tar_cm }
    cm_list \
    | merge_cm \
    | combine(accession_file) \
    | validate_cms \
    | set { gzip_cm }
  emit:
    gzip_cm = gzip_cm
    tar_cm = tar_cm
}
