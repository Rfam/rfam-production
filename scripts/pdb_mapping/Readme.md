# PDB Mapping Pipeline

## Description

The pipeline is run via the pipelines/pdb_mapping.nf script in a different folder,  but all the scripts called by that pipeline are here.

For sequences that have a 3D structure in PDB, we generate a mapping between PDB and Rfam. This pipeline automate this process but the steps are as follows: 
- Get the Rfam.cm file from Rfam CURRENT FTP folder and the latest pdb_seqres file from PDB.
  - The Rfam.cm file needs to compressed with cmpress.
-  Tidy up the Rfam.cm file by removing the "illegal" characters and removing protein sequences. 
-  Run cmscan on all of the PDB sequences. This gives an output ranked list of the CMs of the families with the most significant matches to the sequences. 
-  Import this to the RfamLive database so that you can view the results. 
-  Run the clan competition step. This is a quality assurance measure, run with the aim of reducing redundant hits of families belonging to the same clan.
-  Get the IDs of the families that have been updated with £D information, and the respective PDB IDs.
-  Update the FTP site with the pdb_full_region file. 
-  Update the search index for the website. 
   -  This requires the running of the XML dumper and validator scripts. 
- Update the Release and Web Production databases.
  - We update only the pdb_full_region_table, we do not want to sync these DBs with all of RfamLive outside of release time. 
  - It is necessary to re-run clan competition on the release database. 

This is the end of the PDB Mapping pipeline that uses scripts in the rfam-production repo. We now integrate with the rfam-3d-seed-alignments repo to add the 3D information to the SEED files. On completion of the workflow, a notification is sent to the Rfam slack channel with a summary of the updates. 

## Running

This pipeline will run weekly as a cron job on SLURM. If you would like to run it this can be done by using the batch script like so:

```sbatch scripts/pdb/pdb_mapping.batch```

Or by running the nextflow script:

```nextflow run pipelines/pdb/pdb_mapping.nf```

If you need to run this on an executor other than slurm (e.g. LSF) then update the params in `pipelines/pdb/local.config`

## Notes
The Slack token can be found in `rfam-production/config/rfam_local.py`. If you are testing it may be a good idea to change this to your personal token so as to not flood the Rfam SLack channel with config are up to date. 
The most common issue with the pipeline is that it will not run all the way through past the point of running the 3D seed alignment. If this is the case or will simply need to run the 3D alignment script as per instructions in that repo, as opposed to re-running this whole script. It would be nice if this was improved to have better error handling and/or a retry mechanism. 



