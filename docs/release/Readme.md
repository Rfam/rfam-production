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

### Setup Perl environment

Make sure the `~/.bashrc` file includes the following command:

```
source /path/to/rfam_rh74/rfamrc
```

**Note:** The `rfamrc` file sets up several environment variables, including `PATH`.

### Useful rfamprod aliases to easily move to various locations on the cluster

1. Become `rfamprod` user:

    ```
    become rfamprod
    ```

2. Choose one of the following aliases to move to a specific location on EBI cluster:

    ```
    cd_rel - move to release working directories
    cd_rfamseq - move to Rfamseq location
    cd_rfam - move to rfam-family-pipeline repo
    cd_rh7 - move to Perl libraries
    cd_code - move to rfam-production repo
    cd_main - move to the main Rfam production directory
    ```

---

## Start the release pipeline

1. Update the release version in scripts/release/workflows/nextflow.config 

2. Start the pipeline
```
nextflow run scripts/release/workflows/release_pipeline.nf 
```

This will begin a pipeline that runs the below workflows, in order. If you wish to run these workflows individually, the commands are outlined below. 

## Generate annotated files

This workflow will: 
- Export SEED and CM files
- Generate CM archive zip
- Create tar file 
- Load the SEED and CM files into rfam_live 

```
nextflow run scripts/release/workflows/annotated_files.nf
```

---

## Update PDB mapping

1. Run the PDB mapping pipeline

    The pipeline will update PDB mapping, and then will update the FTP file in nfs/ftp/pub/databases/Rfam/.preview, update the Rfam text search, and begin the process of updating the website database. 

    ```
    nextflow run pdb_mapping/pdb_mapping.nf -profile cluster
    ```

2. Sync the web production databases

    The table pdb_full_region has already been updated in RfamRel database (as part of the pdb_mapping.nf pipeline). This must be synced to the production databses (FB1 and PG). Currently this must be done manually, using the follwoing commands:

    ```
    become mysql-rel-4442
    yes | sync-mysql-fb --dc=FB1
    yes | sync-mysql-fb --dc=PG
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
nextflow run scripts/release/workflows/view_process.nf
```

---

## Clan competition

Clan competition is a quality assurance measure ran as a pre-processing release step aiming to reduce redundant hits of families belonging to the same clan.

```
nextflow run scripts/release/workflows/clan_competition.nf
```

---

## Prepare rfam_live for a new release

This workflow runs the following scripts : 
- [populate_rfamlive_for_release.py](https://github.com/rfam/rfam-production/blob/master/scripts/release/populate_rfamlive_for_release.py)
- [make_rfam_keywords_table.pl](https://github.com/Rfam/rfam-family-pipeline/blob/master/Rfam/Scripts/jiffies/release/make_rfam_keywords_table.pl)
- [updateTaxonomyWebsearch.pl](https://github.com/Rfam/rfam-family-pipeline/blob/master/Rfam/Scripts/jiffies/updateTaxonomyWebsearch.pl) 

```
nextflow run scripts/release/workflows/prepare_rfam_live.nf
```

---

## Generate FTP files

This workflow will generate:
- annotated tree files 
- `Rfam.full_region` file
- `Rfam.pdb` file
- `Rfam.clanin` file
- `fasta_files` folder

```
nextflow run scripts/release/workflows/generate_ftp_files.nf
```

### Generate `database_files` folder
The following steps must be performed manually.

The option `--tab` of `mysqldump` is used to generate dumps in both `.sql` and `.txt` formats, where:
 - `table.sql` includes the MySQL query executed to create a table
 - `table.txt` includes the table contents in tabular format

The main `rfam_live` database is running with the `secure-file-priv` setting so it is not possible to export using the `--tab` option directly. A workaround is to copy dump files on a laptop, and use a local MySQL database for export.

1. Create a new database dump using mysqldump:

    ```
    mysqldump -u <user> -h <host> -P <port> -p --single-transaction --add-locks --lock-tables --add-drop-table --dump-date --comments --allow-keywords --max-allowed-packet=1G --tab=/tmp/database_files rfam_local
    ```

2. Run the [database_file_selector.py](https://github.com/Rfam/rfam-production/blob/master/scripts/release/database_file_selector.py) python script to create the subset of `database_files` available on the FTP:

    ```
    cd /releaseX/database_files
    gzip *.txt
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

---

## Generate `rfam2go` file

```
nextflow run scripts/release/workflows/rfam2go.nf
```

## Stage RfamLive for a new release

This workflow will generate a new MySQL dump to replicate the database on REL and PUBLIC servers. Then the MySQL instance is restored on REL and PUBLIC. 

```
nextflow run scripts/release/workflows/stage_rfam_live.nf
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
nextflow run scripts/release/workflows/update_text_search_dev.nf
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

## Create new json dumps for RNAcentral

### Create a new region export from Rfam using [Rfam2RNAcentral.pl](https://github.com/Rfam/rfam-family-pipeline/blob/master/Rfam/Scripts/export/Rfam2RNAcentral.pl):

1. Extract SEED regions

    ```
    perl Rfam2RNAcentral.pl SEED > /path/to/relX/rnacentral/dir/rfamX_rnac_regions.txt
    ```

2. Extract FULL regions

    ```
    perl Rfam2RNAcentral.pl FULL >> /path/to/relX/rnacentral/dir/rfamX_rnac_regions.txt
    ```

3. Split regions into smaller chunks using basic linux `split` command

    ```
    mkdir chunks &&
    cd chunks &&
    split -n 3000 /path/to/relX/rnacentral/dir/rfamX_rnac_regions.txt rnac_ --additional-suffix='.txt'
    ```

    `-n:` defines the number of chunks to generate (2000 limit on LSF)

    **Notes:**
    - This command will generate 3000 files named like `rnac_zbss.txt`
    - Use `-l 500` option for more efficient chink size

4. Create a copy of the fasta files directory

    ```
    mkdir fasta_files && cd fasta_files &&
    wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/* .
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
    wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.seed.gz && gunzip Rfam.seed.gz
    esl-sfetch --index Rfam.seed
    ```

7. Launch a new json dump using [rnac2json.py](https://github.com/Rfam/rfam-production/blob/master/scripts/export/rnac2json.py):


    - `fasta file directory:` Create a new fasta files directory
    - `json files directory:` Create a new json files output directory


    ```
    python rnac2json.py --input /path/to/chunks --rfam-fasta /path/to/fasta_files --outdir /path/to/json_files
    ```
