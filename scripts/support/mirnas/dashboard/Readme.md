# Dashboard for synchronising miRBase and Rfam

## Update the main Dashboard sheet

The main Dashboard lists the status of all miRBase seed alignments.

1. Get Google Sheet `document_id` and `sheet_id` from the URL:

    `https://docs.google.com/spreadsheets/d/<document_id>/edit#gid=<sheet_id>`

    The document must be shared with `Anyone who has the URL`

1. Update alignments in the `mirbase-seeds` folder, if needed (for example, after a new miRBase release)

1. Update file `microrna_progress.py` with new or updated families (for example, after a new Rfam release)

1. Run `python mirbase_dashboard.py <document_id> <sheet_id>`

    The script will show a URL where the output file can be accessed

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

