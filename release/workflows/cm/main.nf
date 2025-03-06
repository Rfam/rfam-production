process GENERATE_CM_FILE {
  tag "${acc}"
  maxForks 50

  input:
  tuple val(acc), path(seed)

  output:
  path("${acc}.cm")

  """
  writeAnnotatedCM.pl '$acc'
  mv '${acc}.CM' initial.cm
  grep -v DESC initial.cm > '${acc}.nodesc.cm'
  seed-desc-to-cm.pl '$seed' '${acc}.nodesc.cm' > '${acc}.cm'
  """
}

process MERGE_CMS {
  input:
  tuple path('family*.cm'), path(accession_file)

  output:
  path('Rfam.cm'), emit: cm
  path('Rfam.cm.gz'), emit: cm_gz

  """
  find . 'family*.cm' | xargs -I {} cat {} > Rfam.cm
  cmstat Rfam.cm > cmstat.txt
  cm_check.py '${accession_file}' cmstat_file.txt Rfam.cm
  gzip -k Rfam.cm
  """
}

process TAR_CMS {
  input:
  path(cms)

  output:
  path("Rfam.tar.gz")

  """
  tar -czvf Rfam.tar.gz RF0*.cm
  """
}

workflow GENERATE_CM {
  take:
    accessions
    families
    rfam_seed
  emit:
    cm_tar
  main:
    families \
      | combine(rfam_seed) \
      | GENERATE_CM_FILE \
      | collect \
      | set { all_cms }

    all_cms | combine(accessions) | MERGE_CMS | set { cm_gzip }
    all_cms | TAR_CMS | set { cm_tar }

  publish:
    cm_tar >> 'cms'
    cm_gzip >> 'cms'
}
