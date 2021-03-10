# Rfam release process

The Rfam release process involves using tools from the two main GitHub repos

 1. Python repo: `rfam-production`
 2. Perl repo:  `rfam-family-pipeline`

 **Note:** Python code can run locally, but requires EBI VPN connection

## Setup the Python environment

Use virtualenv to install the [requirements](../requirements.txt) locally. Alternatively, run the following on the cluster:

```
become rfamprod
cd_code
source env2/bin/activate
```

## Setup the Perl environment

Modify the `~/.bashrc` file to include the following command:

```
source /path/to/rfam_rh74/rfamrc
```

**Note:** The `rfamrc` file sets up several env variables including `PATH`

## Useful rfamprod aliases to easily move to various locations on the cluster

1. Become `rfamprod` user by executing:

```
become rfamprod
```

2. Choose one of the following aliases to move to a specific location:

```
cd_rel - move to release working directories
cd_rfamseq - move to Rfamseq location
cd_rfam - move to rfam-family-pipeline repo
cd_rh7 - move to Perl libraries
cd_code - move to rfam-production repo
```

**Note:** Requires access to the EBI cluster

---

## Clan competition

Clan competition is a quality assurance measure ran as a pre-processing release step aiming to reduce redundant hits of families belonging to the same clan.

### Preparing for clan competition

1. Create a destination directory for clan competition required files

```
mkdir ~/releaseX/clan_competition
```

2. Generate clan files using [clan_file_generator.py](https://github.com/Rfam/rfam-production/blob/release-14.4/scripts/release/clan_file_generator.py):

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

---

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

---

## Updating PDB sequence file

Updating the `pdb_full_region` table on `rfam_live` depends on the PDBe auto-generated fasta file
and `view` processes.

:information_source: See confluence note

---

## Running view processes

:warning: Requires cluster access and updating the PDB fasta file

1. Create a list of tab separated family accessions and their corresponding uuids using the following query:


```
select rfam_acc, uuid
from _post_process
where status='PEND';
```

2. Update the PDB sequence file and the path in the [PDB plugin](https://github.com/Rfam/rfam-family-pipeline/blob/master/Rfam/Lib/Bio/Rfam/View/Plugin/PDB.pm)

3. Launch view processes on the EBI cluster:

```
python job_dequeuer.py --view-list /path/to/rfam_uuid_pairs.tsv --dest-dir /path/to/destination/directory
```

`--view-list:` A list with tab separated rfam_acc, uuid pairs to run view processes on
`--dest-dir:` The path to the destination directory to generate shell scripts and log to

---

## Launching PDB fasta file scan directly

:information_source: The view processes

1. Scan the PDB fasta file using the newly created `Rfam.cm` on FTP:

```
cmscan -o /path/to/release/dir/PDB_RFAM_X_Y.out --tblout /path/to/release/dir/PDB_RFAM_X_Y.tbl --cut_ga Rfam.cm /path/to/release/dir/PDB_REL_X_Y.fa
```

2. Parse the `PDB_RFAM_X_Y.tbl` and generate a `pdb_full_region` dump (`.txt`) using [infernal_2_pdb_full_region.py](https://github.com/Rfam/rfam-production/blob/master/scripts/processing/infernal_2_pdb_full_region.py):

```
python infernal_2_pdb_full_region.py --tblout /path/to/pdb_full_region.tbl --dest-dir /path/to/dest/dir
```
**NOTE:** This script will generate a new `pdb_full_region_date.txt` file for import to `rfam_live`

3. Truncate `pdb_full_region` table by executing the following query:

```
Truncate table pdb_full_region;
```

4. Import the data in the `.txt` dump to `rfam_live` (see step 2)

```
ADD COMMAND HERE
```

5. Clan compete the hits as described under `Clan competition` section (use PDB option)


### Load SEED and CM files to `rfam_live`:

This enables the SEED and CM download directly from the Rfam website

```
perl populateAnnotatedFiles.pl RFXXXXX ~/path/to/CMs ~/path/to/SEEDs
```

**Note:** This step requires the SEED and CM FTP files

---

## Stage RfamLive for a new release

### Creating MySQL dumps using `mysqldump`

1. Create a new MySQL dump to replicate the database on REL and PUBLIC servers:

```
mysqldump -u username  -h hostname -P port -p --single-transaction --add-locks --lock-tables --add-drop-table --dump-date --comments --allow-keywords --max-allowed-packet=1G rfam_live > rfam_live_relX.sql
```

2. Create a MySQL dump of a single table:

```
mysqldump -u username  -h hostname -P port -p --single-transaction --add-locks --lock-tables --add-drop-table --dump-date --comments --allow-keywords --max-allowed-packet=1G database_name table_name > table_name.sql
```

### Restore a MySQL database instance on a remote server:

1. Move to the directory where the new MySQL dump is located

```
cd /path/to/rfam_live_relX.sql
```

2. Connect to the MySQL server (e.g. PUBLIC)

```
mysql -u username  -h hostname -P port -p
```

3. Create a new MySQL Schema

```
Create schema rfam_X_Y;
```

4. Select schema to restore MySQL dump to

```
Use rfam_X_Y;
```

5. Restore the database in preparation for a new release

```
source rfam_live_relX.sql
```

---

## Generate FTP files

### Generate annotated SEED files:

Export SEED files from SVN using [writeAnnotatedSeed.pl:](https://github.com/Rfam/rfam-family-pipeline/blob/master/Rfam/Scripts/jiffies/writeAnnotatedSeed.pl):

```
perl writeAnnotatedSeed.pl RFXXXXX
```

alternatively use [generate_ftp_files.py](https://github.com/Rfam/rfam-production/blob/release-14.4/scripts/export/generate_ftp_files.py):

```
python generate_ftp_files.py -f /path/to/rfam_accession_list.txt --seed --dest-dir /path/to/seed/files/dest/dir
```

### Generate annotated CM files:

1. Generate a plain CM file using [writeAnnotatedCM.pl:](https://github.com/Rfam/rfam-family-pipeline/blob/master/Rfam/Scripts/jiffies/writeAnnotatedCM.pl):

```
perl writeAnnotatedCM.pl RFXXXXX
```

alternatively use [generate_ftp_files.py](https://github.com/Rfam/rfam-production/blob/release-14.4/scripts/export/generate_ftp_files.py):

```
python generate_ftp_files.py -f /path/to/rfam_accession_list.txt --cm --dest-dir /path/to/CM/files/dest/dir
```

2. Rewrite CM file and descriptions from SEED using [seed-desc-to-cm.pl:](https://github.com/Rfam/rfam-family-pipeline/blob/master/Rfam/Scripts/jiffies/seed-desc-to-cm.pl):

```
perl seed-desc-to-cm.pl <SEED file with DESC> <CM file to add DESC to> > RFXXXXX.cm";
```

### Generate annotated Tree files:

1. Generate new tree files for the release using [writeAnnotatedTree.pl](https://github.com/Rfam/rfam-family-pipeline/blob/master/Rfam/Scripts/jiffies/writeAnnotatedTree.pl):

```
perl writeAnnotatedTree.pl RFXXXXX
```

alternatively use [generate_ftp_files.py](https://github.com/Rfam/rfam-production/blob/release-14.4/scripts/export/generate_ftp_files.py):

```
python generate_ftp_files.py -f /path/to/rfam_accession_list.txt --tree --dest-dir /path/to/tree/files/dest/dir
```

### Export rfam2go

1. Create a new `rfam2go` export by running [rfam2go.pl](https://github.com/Rfam/rfam-family-pipeline/blob/master/Rfam/Scripts/export/rfam2go.pl):

```
rfam2go.pl > /path/to/dest/dir/rfam2go
```

2. Create `md5` checksum of rfam2go file:

```
cd rfam2go &&
md5sum * > md5.txt
```

### Generate `Rfam.full_region` file

To generate the `Rfam.full_region` file execute the following query and dump in a `.txt` file:

```
select rfam_acc, rfamseq_acc, seq_start, seq_end, bit_score, evalue_score, cm_start, cm_end, truncated, type
from full_region
where is_significant=1
order by rfam_acc;
```

### Generate the `Rfam.pdb` file

To generate the `Rfam.pdb` file execute the following query and dump in a `.txt` file:

```
```

### Generate `Rfam.clanin` file

Generate a new `Rfam.clanin` file using [clanin_file_generator.py](https://github.com/Rfam/rfam-production/blob/release-14.4/scripts/release/clanin_file_generator.py):

```
python clanin_file_generator.py --dest-dir /path/to/destination/directory
```

**Note:** Only required if new Clans have been added or existing ones updated

---

## Prepare new data dumps to enable Rfam Text Search

**Requirements:**

The directory for the Text Search index must have the following structure:

-
- families
- clans
- genomes
- motifs
- full_region


1. Move to the main `rfam-production` repo and setup the django environment:

```
source django_settings.sh
```

2. Create an output directory for a new release index
---

## Create new json dumps for import to RNAcentral

1. Create a new region export from Rfam using [Rfam2RNAcentral.pl](https://github.com/Rfam/rfam-family-pipeline/blob/master/Rfam/Scripts/export/Rfam2RNAcentral.pl):

```
perl Rfam2RNAcentral.pl
```

2. Split regions into smaller chunks using basic linux `split` command:

```
ADD COMMAND HERE
```

3. Launch a new json dump using [rnac2json.py]():

```
ADD COMMAND HERE
```

---


---











