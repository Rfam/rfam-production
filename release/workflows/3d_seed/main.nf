process CHECK_3D_ANNOTATION {
  tag { "$acc" }

  input:
  tuple val(acc), path(seed)

  output:
  path("found"), optional: true

  // A family has a 3D seed if it has a GC annotation that ends in _SS, as that is
  // how Rfam marks 3D structures.
  """
  grep '^#=GC' '$seed' | awk '{ print \$3 }' > annotations
  if grep -c '_SS\$' ; then
    cat $seed > found
  fi
  """
}

process MERGE_3D {
  input:
  path("family*.seed")

  output:
  path("Rfam.3d.seed.gz"), emit: seed_gz

  """
  find . -name 'family*.seed' | xargs -I {} cat {} > Rfam.3d.seed
  gzip Rfam.3d.seed
  """
}

workflow GENERATE_3D_SEED {
  take:
    seeds
  main:
    seeds \
    | CHECK_3D_ANNOTATION \
    | collect \
    | MERGE_3D
}
