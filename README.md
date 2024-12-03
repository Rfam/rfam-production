# Rfam Production Code

### About
This repository provides a set of tools related to the Rfam production pipeline.
The collection includes script for data export, database processing, visualisation
and validation, as well as detailed [Rfam release documentation](./release/docs/Release Procedure.md).

### Organisation

This is a mono-repo of several different related tools. Generally these tools
are [nextflow] pipelines.

* **genome-download**: The [nextflow] pipeline for updating genomes in Rfamseq.
* **Export:** Data export tools
* **Processing:** Database processing tools
* **View:** Rfam family view process related scripts
* **Validation:** Data validation scripts
* **Utils:** Support modules/utilities

### Using with Docker

* build container
  ```
  docker build -t rfam-production .
  ```

* open bash interactive shell
  ```
  docker run -v `pwd`:/rfam/rfam-production -it rfam-production bash
  ```

* run a command inside the container
  ```
  docker run -v `pwd`:/rfam/rfam-production -it rfam-production python scripts/export/rfam_xml_dumper.py --out /rfam --type C
  ```

### Installation

```
# install Python dependencies
pip install -r requirements.txt

# set up Python path
export PYTHONPATH=/path/to/the/project

# create and edit configuration file (excluded from version control)
cp config/rfam_local_template.py config/rfam_local.py
```

### Contact us
Feel free to submit an [issue](https://github.com/Rfam/rfam-production/issues) if you have any questions or suggestions.

[nextflow]: <https://www.nextflow.io/>
