# Update and Create microRNA families

## Update miRNA families

Most of these steps use a CSV file as input. This CSV contains the miRBase ID, the Rfam accession number, and the threshold, 
as per the miRBase vs Rfam Dashboard. 
e.g. MIPF0000056__mir-148,RF00248,53
Run the miRNA update pipeline to run the process end-to-end. 
Use the following scripts if you wish to run the processes individually.

1. Check out the family and copy over the SEED:

    `copy_seed.py`

    This takes in a CSV file with the miRNA ID, Rfam accession number, threshold. It runs `rfco.pl`, and copies the SEED to the family directory. 

2. Launch `rfsearch.pl` for the families, given the input CSV:

    `precompute.py`

3. Launch `rfmake.pl` with the threshold per family, given the input CSV:
    `automake.py`

4. Modify the DESC files sequentially:
    `auto_addref.py`

5. Update the fields AU, SE, SS in the DESC files:
    `auto_addref.py`

6. Run the QC checks:
    `auto_rqc.py`

7. Check in the families:
   `auto_rfci.py`

   
## Create new families

1. Create new directory and copy the SEED, for each family:
    `copy_seed.py`
2. Launch `rfsearch.pl` jobs:
    `precompute.py`
3. Launch `rfmake.pl` jobs:
   `automake.py`
4. Generate the correct DESC file, by updating the DESC template:
    `auto_desc_generator.py`
5. Add the microRna references:
    `auto_addref.py`
6. Run the QC checks:
   `auto_rqc.py`
7. Check in the new families, runs `rfnew.pl`:
    `auto_rfnew.py`
