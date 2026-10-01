# PDB Mapping Pipeline

## Description

The pipeline is run via `pipelines/pdb/pdb_mapping.nf`, but all the scripts called by that pipeline are here.

For sequences that have a 3D structure in PDB, we generate a mapping between PDB and Rfam. This pipeline automates the process, and the steps are as follows:
- Get the Rfam.cm file from the Rfam CURRENT FTP folder and the latest pdb_seqres file from PDB.
  - The Rfam.cm file needs to be compressed with cmpress.
- Tidy up the PDB sequences by keeping only nucleic acids and replacing "illegal" characters.
- Run cmscan on all of the PDB sequences. This gives a ranked list of the CMs of the families with the most significant matches to the sequences.
- Import this to the RfamLive database (`pdb_full_region`) so that you can view the results.
- Run the clan competition step. This is a quality assurance measure, run with the aim of reducing redundant hits of families belonging to the same clan.
- Get the IDs of the families that have been updated with 3D information, and the respective PDB IDs (`pdb_families_<date>.txt`).
- Update the FTP site with the pdb_full_region file (`.preview/pdb_full_region.txt.gz`).
- Update the Release and Web Production databases.
  - We update only the pdb_full_region table; we do not want to sync these DBs with all of RfamLive outside of release time.
  - It is necessary to re-run clan competition on the release database.

The website search index is no longer updated here (PR #179); the release pipeline rebuilds it.

On completion, a notification is sent to the Rfam Slack channel with the newest `pdb_families_<date>.txt` report. The message starts with the report date, so a stale date means the current run did not produce a report.

This is the end of the PDB Mapping pipeline that uses scripts in the rfam-production repo. The 3D information is then added to the SEED alignments by `run.sh` in the rfam-3d-seed-alignments repo, which the weekly SLURM batch script runs after this pipeline.

## Running

This pipeline will run weekly as a cron job on SLURM. If you would like to run it this can be done by using the batch script like so:

```sbatch scripts/pdb/pdb_mapping.batch```

The batch script:
- loads a pinned Nextflow version;
- uses `set -eo pipefail`, so the 3D step does not run on a stale mapping if Nextflow fails;
- activates the Python environment, sets `PYTHONPATH` to the repo, and runs the Nextflow pipeline;
- downloads `.preview/pdb_full_region.txt.gz` and runs `run.sh` in rfam-3d-seed-alignments.

Relative `#SBATCH -o`/`-e` paths are resolved against the directory you submit from.

Or by running the nextflow script:

```nextflow run pipelines/pdb/pdb_mapping.nf```

Settings (SLURM executor, paths, memory, cmscan retries) are in `pipelines/pdb/nextflow.config`. Machine-specific overrides can go in a git-ignored `pipelines/pdb/local.config`, passed with `-c`.

## Notes
The Slack token can be found in `rfam-production/config/rfam_local.py`. If you are testing, it may be a good idea to change this to your personal token so as to not flood the Rfam Slack channel with notifications. Alos, so ensure the email addresses in the bacth scripts and the nf config are up to date. 
The most common issue with the pipeline is that it will not run all the way through past the point of running the 3D seed alignment. If this is the case you will simply need to run the 3D alignment script as per instructions in that repo, as opposed to re-running this whole pipeline. It would be nice if this was improved to have better error handling and/or a retry mechanism. 



