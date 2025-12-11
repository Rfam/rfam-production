process CHECK_3D_ANNOTATION {
  tag "$acc"

  input:
  tuple val(acc), path(seed)

  output:
  path "*.seed", optional: true
  //path("found"), optional: true

  // A family has a 3D seed if it has a GC annotation that ends in _SS, as that is
  // how Rfam marks 3D structures.
  script:
  """
  # extract GC annotation column 3
  grep '^#=GC' "$seed" | awk '{ print \$3 }' > annotations

  if grep -q '_SS$' annotations ; then
    #cat $seed > found
    # keep original filename
    cp "$seed" "${acc}.seed"
  fi
  """
}

process MERGE_3D {
  input:
  //path("family*.seed")
  path "*.seed"

  output:
  path("Rfam.3d.seed.gz"), emit: seed_gz

  script:
  """
  #find . -name 'family*.seed' | xargs -I {} cat {} > Rfam.3d.seed
  #gzip Rfam.3d.seed
  if ls *.seed 1> /dev/null 2>&1; then
    cat *.seed > Rfam.3d.seed
  else
    touch Rfam.3d.seed
  fi

  gzip -f Rfam.3d.seed
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
