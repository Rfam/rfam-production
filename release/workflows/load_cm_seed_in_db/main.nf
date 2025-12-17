process PREPARE_MERGED_FILES {
  tag "prepare-merged"

  memory '24 GB'
  time '40m'

  input:
    path seed_gz
    path cm_gz

  output:
    tuple path('Rfam.seed'), path('Rfam.cm')

  script:
  """
  gunzip -c $seed_gz > Rfam.seed
  gunzip -c $cm_gz   > Rfam.cm
  """
}


process LOAD_CM_AND_SEED_ONCE {
  tag "load-cm-seed-db"
  maxForks 1
  memory '64 GB'
  time '20h'

  input:
    tuple path(seed_file), path(cm_file)

  output:
    val("done"), emit: completed

  script:
  """
#!/usr/bin/env python3
import os
import sys
import mysql.connector

def load_cm_seed_in_db(rfam_acc, seed_file, cm_file, cursor, cnx):
    seed_gz_file = f"{rfam_acc}.seed.gz"
    os.system(f"esl-afetch {seed_file} {rfam_acc} | gzip > {seed_gz_file}")
    with open(seed_gz_file, "rb") as f:
        seed_gzip = f.read()
    os.remove(seed_gz_file)

    cm_gz_file = f"{rfam_acc}.cm.gz"
    os.system(f"cmfetch {cm_file} {rfam_acc} | gzip > {cm_gz_file}")
    with open(cm_gz_file, "rb") as f:
        cm_gzip = f.read()
    os.remove(cm_gz_file)

    query = "REPLACE INTO _annotated_file (rfam_acc, seed, cm) VALUES (%s, %s, %s)"
    cursor.execute(query, (rfam_acc, seed_gzip, cm_gzip))
    cnx.commit()

def main(seed_file, cm_file):
    cnx = mysql.connector.connect(
        host="${params.db.host}",
        port=${params.db.port},
        user="${params.db.user}",
        password="${params.db.password}",
        database="${params.db.name}"
    )

    cursor = cnx.cursor()

    # IMPORTANT: this must match the content of Rfam.seed / Rfam.cm
    # Use latin-1 encoding to handle non-UTF-8 characters
    with open(seed_file, encoding='latin-1') as f:
        for line in f:
            if line.startswith('#=GF AC'):
                rfam_acc = line.split()[2]
                print(rfam_acc)
                load_cm_seed_in_db(rfam_acc, seed_file, cm_file, cursor, cnx)

    cursor.close()
    cnx.close()
    print("Done")

if __name__ == "__main__":
    # Use the file names directly from Nextflow
    seed_file = "${seed_file}"
    cm_file = "${cm_file}"
    
    if not os.path.exists(seed_file):
        raise Exception(f"Seed file {seed_file} does not exist")
    
    if not os.path.exists(cm_file):
        raise Exception(f"CM file {cm_file} does not exist")
    
    main(seed_file, cm_file)
  """
}


workflow LOAD_CM_AND_SEED {
  take:
    cm_gzip
    seed_gzip

  main:
    PREPARE_MERGED_FILES(seed_gzip, cm_gzip) \
      | LOAD_CM_AND_SEED_ONCE

  emit:
    completed = LOAD_CM_AND_SEED_ONCE.out.completed
}
