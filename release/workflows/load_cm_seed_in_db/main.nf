process LOAD_CM_AND_SEED {
  input:
  path(family_file)
  path(cm_file)
  path(rfam_seed_file)

  """
  load_cm_seed_in_db.py \
    --host ${params.db.live.host} \
    --port ${params.db.live.port} \
    --user ${params.db.live.user} \
    --password ${params.db.live.password} \
    '${family_file}' '${rfam_seed_file}' '$cm_file'
  """
}
