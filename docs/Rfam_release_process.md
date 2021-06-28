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
cd_main - move to the main Rfam production directory
```

**Note:** Requires access to the EBI cluster

---

## Clan competition

Clan competition is a quality assurance measure ran as a pre-processing release step aiming to reduce redundant hits of families belonging to the same clan.

### Preparing for clan competition

1. Create a destination directory for clan competition required files

    ```
    mkdir ~/releaseX/clan_competition
    mkdir ~/releaseX/clan_competition/sorted    
    ```

2. Generate clan files using [clan_file_generator.py](../scripts/release/clan_file_generator.py):

    ```
    python clan_file_generator.py --dest-dir ~/releaseX/clan_competition --clan-acc CL00001 --cc-type FULL
    ```
    `--cc-type:` Clan competition type FULL/PDB

    `--clan-acc:` Clan accession to compete

    `--all:` Compete all Rfam clans

    `-f:` Only compete clans in file


3. Sort clan files in `dest-dir` based on rfamseq_acc (col2) using linux sort command and store then in the `sorted` directory:

    ```
    sort -k2 -t $'\t' clan_file.txt > ~/releaseX/clan_competition/sorted/clan_file_sorted.txt
    ```

    or for multiple files `cd ~/releaseX/clan_competition` and run:

    ```
    for file in ./CL*; do sort -k2 -t $'\t' ${file:2:7}.txt > sorted/${file:2:7}_s.txt; done
    ```

4. Run clan competition using [clan_competition.py](../scripts/processing/clan_competition.py):

    ```
    python clan_competition.py --input ~/releaseX/clan_competition/sorted --full
    ```

    `--input:` Path to ~/releaseX/clan_competition/sorted

    `--pdb:` Type of hits to compete PDB (pdb_full_region table)

    `--full:` Type of hits to compete FULL (full_region table)

---

## Prepare rfam_live for a new release

1. Truncate `pdb_full_region` table (if launching view processes on all Rfam) or delete the entries for families being updated

2. Populate `rfam_live` tables using [populate_rfamlive_for_release.py](../scripts/release/populate_rfamlive_for_release.py):

    ```
    python populate_rfamlive_for_release.py --all
    ```

3. Make keywords using [make_rfam_keywords_table.pl](https://github.com/Rfam/rfam-family-pipeline/blob/master/Rfam/Scripts/jiffies/release/make_rfam_keywords_table.pl):

    ```
    perl make_rfam_keywords_table.pl
    ```

4. Update `taxonomy_websearch` table using [updateTaxonomyWebsearch.pl](https://github.com/Rfam/rfam-family-pipeline/blob/master/Rfam/Scripts/jiffies/updateTaxonomyWebsearch.pl) (:warning: requires cluster access):

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

## Update PDB mapping

This step requires a finalised `Rfam.cm` file with the latest families, including descriptions (see FTP section for instructions).

1. Get PDB sequences in FASTA format

    ```
    wget ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz
    ```

2. Split into individual sequences

    ```
    perl -pe "s/\>/\/\/\n\>/g" < pdb_seqres.txt > pdb_seqres_sep.txt
    ```

3. Remove protein sequences using [collateSeq.pl](https://github.com/Rfam/rfam-production/blob/master/scripts/release/collateSeq.pl):

    ```
    perl collateSeq.pl pdb_seqres_sep.txt
    mv pdb_seqres_sep.txt.noprot pdb_seqres_sep_noprot.fa
    ```

4. Remove extra text from FASTA descriptor line

    ```
    awk '{print $1}' pdb_seqres_sep_noprot.fa > pdb_trimmed_noprot.fa
    ```

5. Replace any characters that are not recognised by Infernal

    ```
    sed -e '/^[^>]/s/[^ATGCURYMKSWHBVDatgcurymkswhbvd]/N/g' pdb_trimmed_noprot.fa > pdb_trimmed_noillegals.fa
    ```

6. Split into 100 files

    ```
    mkdir files
    seqkit split -O files -p 100 pdb_trimmed_noillegals.fa
    ```

7. Run cmscan in parallel

    ```
    cd files
    cmpress Rfam.cm
    bsub -o part_001.lsf.out -e part_001.lsf.err -M 16000 "cmscan -o part_001.output --tblout part_001.tbl --cut_ga Rfam.cm pdb_trimmed_noillegals.part_001.fa"
    ...
    bsub -o part_100.lsf.out -e part_100.lsf.err -M 16000 "cmscan -o part_100.output --tblout part_100.tbl --cut_ga Rfam.cm pdb_trimmed_noillegals.part_001.fa"
    ```

8. Combine results

    ```
    cat *.tbl | sort | grep -v '#' > PDB_RFAM_X_Y.tbl
    ```

9. Convert Infernal output to `pdb_full_region` table using [infernal_2_pdb_full_region.py](https://github.com/Rfam/rfam-production/blob/master/scripts/processing/infernal_2_pdb_full_region.py):

    ```
    python infernal_2_pdb_full_region.py --tblout /path/to/pdb_full_region.tbl --dest-dir /path/to/dest/dir
    ```    

10. Create a new temporary table in `rfam_live`

    ```
    create table pdb_full_region_temp like pdb_full_region;
    ```

11. Import the data in the `.txt` dump into `rfam_live`

    ```
    mysqlimport -u <user> -h <host> -p -P <port> <database> data_dump.txt
    ```

    Alternatively, use `import` option in Sequel Ace or similar tools

11. Examine the newly imported data and compare with `pdb_full_region`

12. Rename `pdb_full_region` to `pdb_full_region_old` and rename `pdb_full_region_temp` to `pdb_full_region`

13. Clan compete the hits as described under `Clan competition` section using the `--PDB` option

14. List new families with 3D structures

    ```
    # number of families with 3D before
    select count(distinct rfam_acc) from `pdb_full_region_old` where is_significant = 1;

    # number of families with 3D after
    select count(distinct rfam_acc) from `pdb_full_region` where is_significant = 1;    

    # new families with 3D
    select distinct rfam_acc
    from `pdb_full_region`
    where
        is_significant = 1
        and rfam_acc not in
            (select distinct rfam_acc
             from `pdb_full_region_temp`
             where is_significant = 1);
    ```

---

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
export MYSQL_PWD=rfam_live_password
bsub -o mysqldump.out -e mysqldump.err "mysqldump -u <user> -h <hostname> -P <port> --single-transaction --add-locks --lock-tables --add-drop-table --dump-date --comments --allow-keywords --max-allowed-packet=1G rfam_live > rfam_live_relX.sql"
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

Export SEED files from SVN using [generate_ftp_files.py](../scripts/export/generate_ftp_files.py):

```
# to export all files (recommended 16GB RAM to checkout large families)
python generate_ftp_files.py --acc all --seed --dest-dir /path/to/seed/files/dest/dir

# to export accessions listed in a file
python generate_ftp_files.py -f /path/to/rfam_accession_list.txt --seed --dest-dir /path/to/seed/files/dest/dir
```

Alternatively, use [writeAnnotatedSeed.pl](https://github.com/Rfam/rfam-family-pipeline/blob/master/Rfam/Scripts/jiffies/writeAnnotatedSeed.pl):

```
perl writeAnnotatedSeed.pl RFXXXXX
```

### Generate annotated CM files:

1. Export plain CM files using [generate_ftp_files.py](../scripts/export/generate_ftp_files.py):

    ```
    # to export all files (recommended 16GB RAM to checkout large families)
    python generate_ftp_files.py --acc all --cm --dest-dir /path/to/seed/files/dest/dir

    # to export accessions listed in a file
    python generate_ftp_files.py -f /path/to/rfam_accession_list.txt --cm --dest-dir /path/to/CM/files/dest/dir
    ```

    alternatively use [writeAnnotatedCM.pl](https://github.com/Rfam/rfam-family-pipeline/blob/master/Rfam/Scripts/jiffies/writeAnnotatedCM.pl):

    ```
    perl writeAnnotatedCM.pl RFXXXXX
    ```

2. Rewrite CM file and descriptions from SEED using [seed-desc-to-cm.pl](https://github.com/Rfam/rfam-family-pipeline/blob/master/Rfam/Scripts/jiffies/seed-desc-to-cm.pl):

    Required files:

    - `$CM_no_desc` - a CM file to add DESC to (could be 1 CM or all CMs in a single file)
    - `$SEED_with_DESC` - a seed file with DESC lines (could be 1 seed or all seeds in a single file)

    ```
    # filter out DESC lines to avoid duplicates as some CMs already have DESC lines
    grep -v DESC $CM_no_desc > Rfam.nodesc.cm
    perl seed-desc-to-cm.pl $SEED_with_DESC Rfam.nodesc.cm > Rfam.cm

    # check that Rfam.cm contains the correct number of families
    cmstat Rfam.cm | grep -v '#' | wc -l

    # check the number of DESC lines - should be 2 * number of families
    grep DESC Rfam.cm | wc -l

    # generate the final archive
    gzip Rfam.cm
    ```

### Generate annotated Tree files:

1. Generate new tree files for the release using [writeAnnotatedTree.pl](https://github.com/Rfam/rfam-family-pipeline/blob/master/Rfam/Scripts/jiffies/writeAnnotatedTree.pl):

```
perl writeAnnotatedTree.pl RFXXXXX
```

alternatively use [generate_ftp_files.py](../scripts/export/generate_ftp_files.py):

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

1. To generate the `Rfam.pdb` file execute the following query and dump in a `.txt` file:

```
select rfam_acc, pdb_id, chain, pdb_start, pdb_end, bit_score, evalue_score, cm_start, cm_end, hex_colour
from pdb_full_region
where is_significant=1
order by rfam_acc;
```

### Generate `Rfam.clanin` file

Generate a new `Rfam.clanin` file using [clanin_file_generator.py](../scripts/release/clanin_file_generator.py):

```
python clanin_file_generator.py --dest-dir /path/to/destination/directory
```

**Note:** Only required if new Clans have been added or existing ones updated


### Generate FTP database_files

The option `--tab` of `mysqldump` is used to generate dumps in both `.sql` and `.txt` formats, where:
 - `table.sql` includes the MySQL query executed to create a table
 - `table.txt` includes the table contents in tabular format

If executed on the server side, `--tab` can point to a specific location where output will be created. Alternatively, create a local copy of the database and
dump the files to `tmp` directory.

1. Create a new database dump using mysqldump:

```
mysqldump -u <user> -h localhost -P <port> -p --single-transaction --add-locks --lock-tables --add-drop-table --dump-date --comments --allow-keywords --max-allowed-packet=1G --tab=/tmp rfam_live_local
```

2. Run the [database_file_selector.py](https://github.com/Rfam/rfam-production/blob/master/scripts/release/database_file_selector.py) python script to create the subset of `database_files` available on the FTP:

```
python database_file_selector.py --source-dir /tmp --dest-dir /releaseX/database_files
```

3. Zip `database_files` directory and copy to remote server:

```
tar -czvf database_files.tar.gz /releaseX/database_files
scp database_files.tar.gz username@remote.host:/some/location
```

4. Restore `database_files` on FTP

```
tar -xzvf /some/location/database_files.tar.gz .
```

### Generate FTP fasta_files

1. Export new fasta files for all families in Rfam using [fasta_file_generator.py](https://github.com/Rfam/rfam-production/tree/master/scripts/export/fasta_file_generator.py):

```
# Submit a new interactive LSF job
bsub -M 10000 -Is $SHELL
# launch the python script to export fasta files
python fasta_file_generator.py --seq-db /path/to/rfamseq.fa --rfam-seed /path/to/releaseX/Rfam.seed --all --outdir /path/to/ftp/fasta_files
```

ALternatively, launch individual LSF jobs using X:

```
ADD COMMAND HERE
```

2. Copy the newly generated `fasta_files` to FTP

---

## Prepare new data dumps to enable Rfam Text Search

### Generate new data dumps

**Requirements:**

The directory for the Text Search index must have the following structure:

- release_note.txt
- families
- clans
- genomes
- motifs
- full_region


1. Move to the main `rfam-production` repo and setup the django environment:

```
source django_settings.sh
```

2. Create an output directory for a new release index and all required subdirectories:

```
mkdir -p relX_text_search/families
```

3. Create a new families XML dump:

```
python rfam_xml_dumper.py --type F --out /path/to/relX_text_search/families
```

4. Validate xml dumps:

```
python xml_validator.py --input /path/to/relX_text_search/families --log
```

`--input`: Single XML dump or a directory
 `--log`: Generates a log file with all XMl dumps failing validation


Follow the same process to generate  XML dumps for all data types:
 - families
 - clans
 - motifs
 - genomes
 - full_region

For more information on launching XML dumps as LSF jobs see [lsf_rfam_xml_dumper.sh](../scripts/export/lsf_rfam_xml_dumper.sh)

:information_source: More detailed information available on confluence

### Index data on **dev** and **prod**

1. Change directory to `text_search` on the cluster:

```
cd_main && cd search_dumps
```

2. Index data on `dev`:

```
unlink rfam_dev
ln -s /path/to/xml/data/dumps rfam_dev
```

3. Index data on `prod`:

```
unlink current_release
ln -s /path/to/xml/data/dumps current_release
```
---

## Create new json dumps for import to RNAcentral


### Create a new region export from Rfam using [Rfam2RNAcentral.pl](https://github.com/Rfam/rfam-family-pipeline/blob/master/Rfam/Scripts/export/Rfam2RNAcentral.pl):

1. Extract SEED regions:
```
perl Rfam2RNAcentral.pl SEED > /path/to/relX/rnacentral/dir/rfamX_rnac_regions.txt
```

2. Extract FULL regions:
```
perl Rfam2RNAcentral.pl FULL >> /path/to/relX/rnacentral/dir/rfamX_rnac_regions.txt
```

3. Split regions into smaller chunks using basic linux `split` command:

```
mkdir chunks &&
cd chunks &&
split -n 3000 /path/to/relX/rnacentral/dir/rfamX_rnac_regions.txt rnac_ --additional-suffix='.txt'
```

`-n:` defines the number of chunks to generate (2000 limit on LSF)

**Notes:**
- This command will generate 3000 files named like `rnac_zbss.txt`
- Use `-l 500` option for more efficient chink size

4. Create a copy of the fasta files directory:

```
mkdir fasta_files && cd fasta_files &&
wget http://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/* .
```
**NOTE:** Rfam2RNAcentral regions need to match the exact same release the `fasta_files` were created for.
          Use the Public MySQL database if `rfam-live` have been updated

5. Unzip all fasta files and index using `esl-sfetch`:

```
gunzip * .gz
for file in ./RF*.fa; do esl-sfetch --index $file; done
```

6. Copy the `Rfam.seed file from the FTP to the fasta_files directory and index using `esl-sfetch`:

```
wget http://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.seed.gz && gunzip Rfam.seed.gz
esl-sfetch --index Rfam.seed
```

7. Launch a new json dump using [rnac2json.py](https://github.com/Rfam/rfam-production/blob/master/scripts/export/rnac2json.py):


- `fasta file directory:` Create a new fasta files directory
- `json files directory:` Create a new json files output directory


```
python rnac2json.py --input /path/to/chunks --rfam-fasta /path/to/fasta_files --outdir /path/to/json_files
```

---
