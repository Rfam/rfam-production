process GENERATE_CM_FILE {
  tag "${acc}"
  maxForks 50
  memory { acc == "RF00005" ? 10.GB 2.GB }

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
  find . -name 'family*.cm' | xargs -I {} cat {} > Rfam.cm
  gzip Rfam.cm
  """
}

workflow GENERATE_CM {
  take:
    seed_alignments
  emit:
    cm_gzip
    all_cms
  main:
    seed_alignments | GENERATE_CM_FILE | set { all_cms }

    all_cms | collect | MERGE_CMS | set { cm_gzip }
}
