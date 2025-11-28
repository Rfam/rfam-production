process BUILD {
  input:
  path(query)

  output:
  path("Rfam.clanin")

  script:
  """
  set -euo pipefail
  
  tmp=\$(mktemp)
  cleanup() { rm -f "\$tmp"; }
  trap cleanup EXIT

  mysql -s \
    --host=${params.db.host} \
    --port=${params.db.port} \
    --user=${params.db.user} \
    --database=${params.db.name} \
    --password=${params.db.password} \
    < "${query}" > "\$tmp"

  # verify output is not empty
  if [ ! -s "\$tmp" ]; then
    echo "ERROR: Query returned no results" >&2
    exit 1
  fi
  
  mv "\$tmp" Rfam.clanin
  trap - EXIT
  """
}

workflow GENERATE_CLANIN {
  emit:
  clanin
  main:
  channel.fromPath("${moduleDir}/sql/clanin.sql") | BUILD | set { clanin }
}
