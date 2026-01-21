process BUILD {
  tag "${acc}"
  maxForks 50
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 3
  memory { 
    [2.GB, 30.GB, 120.GB][task.attempt - 1] ?: 120.GB
  }

  input:
  tuple val(acc), path(seed)

  output:
  path("${acc}.cm")

  script:
  """
  writeAnnotatedCM.pl '${acc}'
  mv '${acc}.CM' initial.cm
  grep -v DESC initial.cm > '${acc}.nodesc.cm'
  seed-desc-to-cm.pl '${seed}' '${acc}.nodesc.cm' > '${acc}.cm'
  """
}

process MERGE_CMS {
  input:
  path('family*.cm')

  output:
  path('Rfam.cm.gz'), emit: cm_gz

  script:
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
  seed_alignments | BUILD | set { all_cms }
  all_cms | collect | MERGE_CMS | set { cm_gzip }
}
