# Update and Create microRNA families

## Update miRNA families

Most of these steps use a TSV file as input. This TSV contains the miRBase ID, the Rfam accession number, and the threshold, 
as per the miRBase vs Rfam Dashboard. 
e.g. MIPF0000056__mir-148,RF00248,53

Use the following scripts if you wish to run the processes individually.

Ensure the directories in `mirna_config.py` are as expected.

1. Check out the family and copy over the SEED:

    `copy_seed.py --input <file>`

    This takes in a TSV file with the miRNA ID, Rfam accession number, threshold. It runs `rfco.pl`, and copies the SEED to the family directory. 

2. Launch `rfsearch.pl` for the families, given the input TSV:
    `auto_rfsearch.py --input <file>`

3. Launch `rfmake.pl` with the threshold per family, given the input TSV:
    `auto_rfmake.py --input <file>`

4. Modify the DESC files sequentially:
    `auto_addref.py --input <file> --sequential`

5. Update the fields AU, SE, SS in the DESC files:
    `update_desc.py --input <file>`

6. Run the QC checks:
    `auto_rqc.py --input <file>`

7. Check in the families:
   `auto_rfci.py --input <file>`

   
## Create new miRNA families

1. Create a new directory and copy the SEED, for each family:
    This takes in a TSV file with the miRNA ID, threshold. It runs `rfco.pl`, and copies the SEED to the family directory. 
    Use argument `--new` to specify that these are new miRNA families. 
    
    `copy_seed.py --new --input <file>`

2. Launch `rfsearch.pl` jobs:
    Use argument `--new` to specify that these are new miRNA families.

    `auto_rfsearch.py --new --input <file>`

3. Launch `rfmake.pl` jobs:
   
4. `auto_rfmake.py --new --input <file>`

5. Generate the correct DESC file, by updating the DESC template:
    Only needs to be run for new families, families to be updated will already have a DESC file

    `auto_desc_generator.py --input <file>`

6. Add the microRNA references:
    
    `auto_addref.py --new --input <file> --sequential`

7. Update the fields AU, SE, SS in the DESC files:
    
    `update_desc.py --new --input <file>`

8. Run the QC checks:
   
    `auto_rqc.py --new --input <file>`

9. Check in the new families, runs `rfnew.pl`:
     
    `auto_rfnew.py --input <file>`
