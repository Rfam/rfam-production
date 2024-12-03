# Rfam release process

## Setup

The Rfam release process uses two main GitHub repos:

 1. Python repo: [rfam-production](https://github.org/rfam/rfam-production)
 2. Perl repo:  [rfam-family-pipeline](https://github.org/rfam/rfam-family-pipeline)

 **Note:** Most Python code can run locally, but requires an EBI VPN connection. All Perl code currently runs only on the EBI cluster.

### Setup Python environment

Run the following on the EBI cluster:

```
become rfamprod
cd_code && cd rfam-production
source env/bin/activate
```

Alternatively, use virtualenv to install the [requirements](../requirements.txt) locally.

### Useful rfamprod aliases to easily move to various locations on the cluster

1. Become `rfamprod` user:

    ```
    become rfamprod
    ```

2. Choose one of the following aliases to move to a specific location on EBI cluster:

    ```
    cd_code - move to rfamprd code folder
    cd_rel - move to release working directories
    ```

---

## Start the release process

1. Update the release version in pipelines/release/local.config
2. Create a release_x folder in the release directory

3. Start the pipeline
```
nextflow run scripts/release/workflows/release_pipeline.nf
```

This will begin a pipeline that runs the below workflows, in order. If you wish to run these workflows individually, the commands are outlined below.

**Note:** This is not working correctly, so it is best to run the pipelines individually.

## Generate annotated files

This workflow will:
- Export SEED and CM files
- Generate CM archive zip
- Create tar file
- Load the SEED and CM files into rfam_live

```
nextflow run pipelines/release/annotated_files.nf
```

---

## Update PDB mapping

Run the PDB mapping pipeline
- The pipeline will update PDB mapping, and then will update the FTP file in nfs/ftp/public/databases/Rfam/.preview, update the Rfam text search, and begin the process of updating the website database.
- This script needs some modification to run for release. It usually runs weekly, using the Rfam.cm file from the FTP site. For release, we need it to use the newly created CM file.

    ```
    nextflow run pdb_mapping/pdb_mapping.nf -profile cluster
    ```

<details>
  <summary>Legacy steps for manually updating PDB mapping</summary>
Please note these steps for updating PDB Mapping have been replaced with the introduction of the above PDB mapping pipeline.


This step requires a finalised `Rfam.cm` file with the latest families, including descriptions (see FTP section for instructions).

1. Get PDB sequences in FASTA format

    ```
    wget ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz
    gunzip pdb_seqres.txt.gz
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

    To check that all commands completed, check that every log file contains `Successfully completed`:

    ```
    cat part_*lsf.out | grep Success | wc -l
    ```

8. Combine results

    ```
    cat *.tbl | sort | grep -v '#' > PDB_RFAM_X_Y.tbl
    ```

9. Convert Infernal output to `pdb_full_region` table using [infernal_2_pdb_full_region.py](https://github.com/Rfam/rfam-production/blob/master/scripts/processing/infernal_2_pdb_full_region.py):

    ```
    python infernal_2_pdb_full_region.py --tblout /path/to/pdb_full_region.tbl --dest-dir /path/to/dest/dir
    ```

    The script will generate a file like `pdb_full_region_YYYY-MM-DD.txt`.

10. Create a new temporary table in `rfam_live`

    ```
    create table pdb_full_region_temp like pdb_full_region;
    ```

11. Manually import the data in the `.txt` dump into `rfam_live` using Sequel Ace or similar

    :warning: `mysqlimport` or `LOAD DATA INFILE` do not work with `rfam_live` because of `secure-file-priv` MySQL setting.

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
</details>

---

## Running view processes

```
nextflow run pipelines/release/view_process.nf
```

---

## Clan competition

Clan competition is a quality assurance measure ran as a pre-processing release step aiming to reduce redundant hits of families belonging to the same clan.

```
nextflow run pipelines/release/clan_competition.nf
```

---

## Prepare rfam_live for a new release

This workflow runs the following scripts :
- [populate_rfamlive_for_release.py](https://github.com/rfam/rfam-production/blob/master/scripts/release/populate_rfamlive_for_release.py)
- [make_rfam_keywords_table.pl](https://github.com/Rfam/rfam-family-pipeline/blob/master/Rfam/Scripts/jiffies/release/make_rfam_keywords_table.pl)
- [updateTaxonomyWebsearch.pl](https://github.com/Rfam/rfam-family-pipeline/blob/master/Rfam/Scripts/jiffies/updateTaxonomyWebsearch.pl)

```
nextflow run pipelines/release/prepare_rfam_live.nf
```

---

## Generate FTP files

This workflow will generate:
- annotated tree files
- `Rfam.full_region` file
- `Rfam.pdb` file
- `Rfam.clanin` file
- `fasta_files` folder
- `Rfam.3d.seed.gz`

```
nextflow run pipelines/release/generate_ftp_files.nf
```

### Generate `database_files` folder
The following steps must be performed manually.

The option `--tab` of `mysqldump` is used to generate dumps in both `.sql` and `.txt` formats, where:
 - `table.sql` includes the MySQL query executed to create a table
 - `table.txt` includes the table contents in tabular format

The main `rfam_live` database is running with the `secure-file-priv` setting, so it is not possible to export using the `--tab` option directly. A workaround is to copy dump files on a laptop, and use a local MySQL database for export.

1. Export a the .sql files of the rfam_live database (using Sequel Ace or via command line), then recreate a local copy, for example
   ```
   mysql > use rfam_live_14_8;
   mysql > source /Users/user/Desktop/sql_queries/14_8_copy.sql;
   ```
2. Create a new database dump of the local copy of rfam_love using mysqldump:

    ```
    mysqldump -u <user> -h <host> -P <port> -p --single-transaction --add-locks --lock-tables --add-drop-table --dump-date --comments --allow-keywords --max-allowed-packet=1G --tab=/tmp/database_files rfam_local
    ```

3. Run the [database_file_selector.py](https://github.com/Rfam/rfam-production/blob/master/scripts/release/database_file_selector.py) python script to create the subset of `database_files` available on the FTP:

    ```
    cd /tmp/database_files
    gzip *.txt
    python database_file_selector.py --source-dir /tmp/database_files --dest-dir /releaseX/database_files
    ```

4. Zip `database_files` directory and copy to remote server:

    ```
    tar -czvf database_files.tar.gz /releaseX/database_files
    scp database_files.tar.gz username@remote.host:/some/location
    ```

5. Restore `database_files` on FTP

    ```
    tar -xzvf /some/location/database_files.tar.gz /path/to/release/ftp/database_files
    ```

---

## Generate `rfam2go` file

```
nextflow run pipelines/release/rfam2go.nf
```

## Stage RfamLive for a new release

This workflow will generate a new MySQL dump to replicate the database on REL and PUBLIC servers. Then the MySQL instance is restored on REL and PUBLIC.

```
nextflow run pipelines/release/stage_rfam_live.nf
```

---

## Update Rfam text search

Description of indexed fields: https://www.ebi.ac.uk/ebisearch/metadata.ebi?db=rfam

### Generate new data dumps

**Requirements:**

The directory for the text search index must have the following structure:

- release_note.txt
- families
- clans
- genomes
- motifs
- full_region

This workflow will also index the data on dev. Once this is successful, it is necessary to index the data on **prod**

```
nextflow run pipelines/release/update_text_search_dev.nf
```

### Index data on **prod**

1. Change directory to `text_search` on the cluster:

    ```
    cd_main && cd search_dumps
    ```

2. Index data on `prod`:

    ```
    unlink current_release
    ln -s /path/to/xml/data/dumps current_release
    ```

The files in `rfam_dev` and `current_release` folders are automatically indexed every night.

---
