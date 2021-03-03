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

1. Create a destination directory for clan competition required files

```
mkdir ~/releaseX/clan_competition
```


2. Generate clan files using [clan_file_generator.py](https://github.com/Rfam/rfam-production/blob/release-14.4/scripts/release/clan_file_generator.py)

```
python clan_file_generator.py --dest-dir ~/releaseX/clan_competition --clan-acc CL00001 --cc-type FULL
```
`--cc-type:` Clan competition type FULL/PDB

`--clan-acc:` Clan accession to compete

`--all:` Compete all Rfam clans

`-f:` Only compete clans in file


3. Sort clan files in `dest-dir` based on rfamseq_acc (col2) using linux sort command and store then in the sorted directory:

```
sort -k2 -t $'\t' clan_file.txt > ~/releaseX/clan_competition/sorted/clan_file_sorted.txt
```

or for multiple files `cd ~/releaseX/clan_competition` and run:

```
for file in ./CL*; do sort -k2 -t $'\t' ${file:2:7}.txt > sorted/${file:2:7}_s.txt; done
```

4. Run clan competition using [clan_competition.py](https://github.com/Rfam/rfam-production/blob/release-14.4/scripts/processing/clan_competition.py):

```
python clan_competition.py --input ~/releaseX/clan_competition/sorted --full
```

`--input:` Path to ~/releaseX/clan_competition/sorted

`--pdb:` Type of hits to compete PDB (pdb_full_region table)

`--full:` Type of hits to compete FULL (full_region table)

## Prepare rfam_live for a new release
1. Truncate `family_ncbi` table

2. Truncate `pdb_full_region` table (if launching view processes on all Rfam) or delete the entries for families being updated
3. Populate `rfam_live` tables using [populate_rfamlive_for_release.py](https://github.com/Rfam/rfam-production/blob/release-14.4/scripts/release/populate_rfamlive_for_release.py):

```
python populate_rfamlive_for_release.py --all
```
4. Make useful keywords by running:

```
make_rfam_keywords_table.pl
```

5. Update `taxonomy_websearch` table using [updateTaxonomyWebsearch.pl](https://github.com/Rfam/rfam-family-pipeline/blob/master/Rfam/Scripts/jiffies/updateTaxonomyWebsearch.pl) (:warning: requires cluster access):

```
perl updateTaxonomyWebsearch.pl
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

2. Launch view processes on the EBI cluster:

```
python job_dequeuer.py --view-list /path/to/rfam_uuid_pairs.tsv --dest-dir /path/to/destination/directory
```

`--view-list:` A list with tab separated rfam_acc, uuid pairs to run view processes on
`--dest-dir:` The path to the destination directory to generate shell scripts and log to


## Stage RfamLive for a new release

Create a new MySQL dump to replicate the database on REL and PUBLIC servers

```
mysqldump -u username  -h hostname -P port -p --single-transaction --add-locks --lock-tables --add-drop-table --dump-date --comments --allow-keywords --max-allowed-packet=1G rfam_live > rfam_live_relX.sql
```

## Generate FTP files

### Generate SEED files using [writeAnnotatedSeed.pl:](https://github.com/Rfam/rfam-family-pipeline/blob/master/Rfam/Scripts/jiffies/writeAnnotatedSeed.pl)

```
perl writeAnnotatedSeed.pl RFXXXXX
```

### Generate annotated CM files:

1. Generate a plain CM file using [writeAnnotatedCM.pl:](https://github.com/Rfam/rfam-family-pipeline/blob/master/Rfam/Scripts/jiffies/writeAnnotatedCM.pl)

```
perl writeAnnotatedCM.pl RFXXXXX
```

2. Rewrite CM file and descriptions from SEED using [seed-desc-to-cm.pl:](https://github.com/Rfam/rfam-family-pipeline/blob/master/Rfam/Scripts/jiffies/seed-desc-to-cm.pl)

```
perl seed-desc-to-cm.pl <SEED file with DESC> <CM file to add DESC to> > RFXXXXX.cm";
```

### Generate annotated Tree files:

1. Generate new tree files for the release using [writeAnnotatedTree.pl](https://github.com/Rfam/rfam-family-pipeline/blob/master/Rfam/Scripts/jiffies/writeAnnotatedTree.pl):

```
perl writeAnnotatedTree.pl RFXXXXX
```

### Load SEED and CM files to `rfam_live`:

This enables the SEED and CM download directly from the Rfam website

```

```




