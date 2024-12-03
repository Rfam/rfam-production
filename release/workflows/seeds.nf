#!/usr/bin/env nextflow

process family_seed {
  tag "${accession}"
  memory { acc == "RF00005" ? 20.GB : 5.GB }
  // Seems that writeAnnotatedSeed has creates invalid SEED files, missing
  // metadata, if too many process run at once. 2 seems to work correctly, 10
  // does not. Let's go with the lower value
  maxForks 2

  input:
  val(accession)

  output:
  path("${accession}.seed")

  """
  writeAnnotatedSeed.pl $accession
  """
}

process merge_seed {
  publishDir "${params.release_ftp}", mode: "copy"

  input:
  path(seed)

  output:
  path('Rfam.seed.gz')

  """
  cat $seed > Rfam.seed
  esl-afetch --index Rfam.seed
  gzip -c Rfam.seed > Rfam.seed.gz
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

process build_3d_seed {
  publishDir "${params.release_ftp}", mode: "copy"

  input:
  tuple path('Rfam.seed.gz'), path('familes.txt')

  output:
  path('Rfam.3d.seed.gz')

  """
  zcat Rfam.seed.gz > Rfam.seed
  esl-afetch -f Rfam.seed families.txt > Rfam.3d.seed
  gzip Rfam.3d.seed
  """
}

workflow seeds {
  take:
    families
    fam_3d
  main:
    families | family_seed | set { all_seeds }

    all_seeds | collect | merge_seed | set { rfam_seed }
    rfam_seed | combine(fam_3d) | build_3d_seed | set { seed_3d }
  emit:
    rfam_seed = rfam_seed
    seed_3d = seed_3d
}
