// process FETCH_FAMILIES {
//   output:
//   path('families')
//
//   """
//   mysql -s \
//     --host=${params.db.host} \
//     --port=${params.db.port} \
//     --user=${params.db.user} \
//     --database=${params.db.name} \
//     --password=${params.db.password} \
//     <<< "select rfam_acc from family" > families
//   """
// }


process FETCH_FAMILIES {
  output:
  path('families')
  
  script:
  """
  set -euo pipefail
  tmp=\$(mktemp)
  cleanup() { rm -f "\$tmp"; }
  trap cleanup EXIT
  mysql -s \\
    --host=${params.db.host} \\
    --port=${params.db.port} \\
    --user=${params.db.user} \\
    --database=${params.db.name} \\
    --password=${params.db.password} \\
    -e "SELECT rfam_acc FROM family" > "\$tmp"
  
  # verify output is not empty
  if [ ! -s "\$tmp" ]; then
    echo "ERROR: Query returned no results" >&2
    exit 1
  fi
  mv "\$tmp" families
  trap - EXIT
  """
}
