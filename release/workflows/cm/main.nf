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
  path('family*.cm')

  output:
  path('Rfam.cm.gz'), emit: cm_gz

  """
  find . 'family*.cm' | xargs -I {} cat {} > Rfam.cm
  gzip Rfam.cm
  """
}

workflow GENERATE_CM {
  take:
    seed_alignments
  emit:
    cm_gzip
  main:
    seed_alignments | GENERATE_CM_FILE | set { all_cms }

    all_cms | collect | MERGE_CMS | set { cm_gzip }
  publish:
    all_cms >> 'all_cms'
    cm_gzip >> 'cms'
}
