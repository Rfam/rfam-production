# Installation
Create a python environment and install the dependencies defined in cwl-requirements.txt
Docker / Singularity is required to run this workflow!

# Executing the CWL workflow
To run the cwl workflow with a defined input;
Create a file test.yml using the following template:
```yml
sequences:
  class: File
  path: ./RF02567.fa
  format: http://edamontology.org/format_1929
covariance_model_database:
  class: File
  path: ../db/RF02567.cm
```
Then run:
`cwltool --singularity --out <output_directory> nc_rna_workflow.cwl test.yml`

## Executing from python script: 
Run wrapper_scripts/process_dir.py ]

## Custom docker container
A custom docker container can be built using the CWL-Dockerfile.
This will bundle all the rfam-production code into a simple python container,
which is called by many tools in cwl/tools.  
