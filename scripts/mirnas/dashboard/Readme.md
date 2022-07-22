# Dashboard for synchronising miRBase and Rfam

## Update the main Dashboard sheet

The main Dashboard lists the status of all miRBase seed alignments.

1. Get Google Sheet `document_id` and `sheet_id` from the URL:

    `https://docs.google.com/spreadsheets/d/<document_id>/edit#gid=<sheet_id>`

    The document must be shared with `Anyone who has the URL`

1. Update alignments in the `mirbase-seeds` folder, if needed (for example, after a new miRBase release)

1. Update file `microrna_progress.py` with new or updated families (for example, after a new Rfam release)

1. Run `python mirbase_dashboard.py <document_id> <sheet_id>`

    - Use `--nocache` to recompute all overlaps (they are cached by default)
    - The script will show a URL where the output file can be accessed

1. Download the newly generated file to your computer

1. In Google Sheets, create a new sheet by duplicating the current Dashboard

1. Go to File > Import > Upload and select the new file

1. Select `Replace current sheet` and make sure `Convert text to numbers, dates and formulas` is selected

1. Sort the data

    1. Select all columns with data except for columns A and B
    1. Select Data > Sort range > Advanced range sorting options
    1. Click `Data has header row`
    1. Select `Action` or any other column to sort by in the dropdown

1. Optional: copy Conditional formatting across sheets

    1. Select the entire column with Conditional formatting and do Edit > Copy
    1. Select a target column and do Edit > Paste special > Conditional format only

1. Check that the summary columns make sense compared to the previous version

1. Rename old Dashboard sheet and hide it, rename the new Dashboard sheet

## Update the Family updates sheet

1. Get Google Sheet `document_id` and `sheet_id` from the URL as above.

    - The main Dashboard is used to get a list of  families that require updates.

1. Run `python find_family_overlaps.py <document_id> <sheet_id>`

    - The script will show a URL where the output file can be accessed

1. Download the newly generated file to your computer

1. In Google Sheets, create a new sheet by duplicating the current `Family updates` sheet

1. Go to File > Import > Upload and select the new file

1. Select `Replace current sheet` and make sure `Convert text to numbers, dates and formulas` is selected

1. Check that the summary columns make sense compared to the previous version

1. Rename old Dashboard sheet and hide it, rename the new Dashboard sheet

## Generate an archive for 1_SEED families

To facilitate manual review of entries with a single seed sequence, an archive
with relevant text files can be generated as follows:

1. Get Google Sheet `document_id` and `sheet_id` from the URL as above.

1. Run `python process_one_seed_families.py <document_id> <sheet_id`

    - The script will show a URL where the output file can be accessed

## Update `Rfam miRNA without matches` and  `Rfam non-miRNA families matching miRBase` sheets

This step requires a mapping between all sequences from miRBase and the current Rfam models.
The mapping should be recomputed after a new Rfam release or when a new set of covariance models is available.
A precomputed file `mirbase-cmscan.tblout` is included in the repository.

1. (Optional) Update file `mirbase-cmscan.tblout` as follows (takes ~2.5 hours on a laptop):

    ```
    wget https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/sequences/by-database/mirbase.fasta .
    wget https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz .
    gunzip Rfam.cm.gz && cmpress Rfam.cm
    cmscan -o cmscan.output.txt --tblout mirbase-cmscan.tblout --cut_ga --rfam Rfam.cm mirbase.fasta
    ```

    - Commit new file for the record

2. Run `python compare_rfam_mirbase.py` to analyse a precomputed file or `python compare_rfam_mirbase.py <tblout>` to analyse a different file

    - The script will print a list of Rfam microRNA families without miRBase hits. At the time of writing it includes only 2 families.
    - The script will show a URL where the output file with Rfam non-microRNA families matching miRBase can be accessed

3. Download the newly generated file to your computer

4. In Google Sheets, create a new sheet by duplicating the current `Rfam non-miRNA families matching miRBase` sheet

5. Go to File > Import > Upload and select the new file

6. Select `Replace current sheet` and make sure `Convert text to numbers, dates and formulas` is selected