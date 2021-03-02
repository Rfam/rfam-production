# Rfam release process

## Setup the environment

Use virtualenv to install the [requirements](../requirements.txt) locally. Alternatively, run the following on the cluster:

```
become rfamprod
cd_code
source env2/bin/activate
```

## Clan competition

Clan competition is a quality assurance measure ran as a pre-processing release step aiming to reduce redundant hits of families belonging to the same clan. 

### Preparing for clan competition

1. Generate clan files using [clan_file_generator.py](https://github.com/Rfam/rfam-production/blob/release-14.4/scripts/release/clan_file_generator.py)

```
python clan_file_generator.py --dest-dir /path/to/dest/directory --clan-acc CL00001 --cc-type FULL
```
`--dest-dir:` Destination directory to generate output to

`--cc-type:` Clan competition type **[FULL|PDB]**

`--clan-acc:` Clan accession to compete

`--all:` Compete all Rfam clans

`-f:` Only compete clans in file

2. Sort clan files in dest-dir based on rfamseq_acc (col2) using linux sort command:

```
sort -k2 -t $'\t\' clan_file.txt > clan_file_sorted.txt
```

3. Run clan competition using [clan_competition.py](https://github.com/Rfam/rfam-production/blob/release-14.4/scripts/processing/clan_competition.py):

```
clan_competition.py 
```

## Prepare rfam_live for a new release
1. Truncate `family_ncbi` table

2. Truncate `pdb_full_region` table (if launching view processes on all Rfam) or delete the entries for families being updated
3. Populate `rfam_live` tables using [populate_rfamlive_for_release.py](https://github.com/Rfam/rfam-production/blob/release-14.4/scripts/release/populate_rfamlive_for_release.py):

```
python populate_rfamlive_for_release.py --all
```
4. Make useful keywords by running (requires cluster access):

```
perl make_rfam_keywords_table.pl
```

## Updating PDB sequence file

:information_source: See confluence notes

## Running view processes

:warning: Requires cluster access

1. Create a list of tab separated family accessions and their corresponding uuids using the following query:

```
select rfam_acc, uuid 
from _post_process 
where status='PEND';
```

2. Launch view processes on the cluster:

```
job_dequeuer.py 
```




---
**NOTE**

---


