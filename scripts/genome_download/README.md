Rfamseq Update
==============

This is a [nextflow](https://www.nextflow.io/) pipeline meant to updated Rfamseq by fetching all required genomes and produce the metadata required for import into the Rfam mysql database.

What is Rfamseq?
----------------

Please read over the docs about [rfamseq](https://docs.rfam.org/en/latest/glossary.html#rfamseq) along with the cited papers for full details.
In short, this is a collection of genomes that Rfam uses to build the 'full hits' and is searched by `rfsearch`.
This set of genomes defines the range of genomes Rfam can find hits to and should be as comprehensive as is reasonable.
Currently, we use [UniProt](https://www.uniprot.org/) reference proteomes to build the genomes.
This is useful because this is a proteome set which aims to solve a similar problem, allows our users to see protein and ncRNA annotations in the same genome and is handled by other people making it easier for us.

This pipeline will produce the files needed for Rfamseq. They are:

- `Rfam_{version}.fa` - The complete database.
- `Rfam_{version}.fa.ssi` - The index for the complete database.
- `rfamseq_{version}_*.fa` - A series of chunks of the fasta above.
- `rev-rfamseq_*.fa` - A series of chunks of the reversed database.

Note that the Rfam database also contains metadata about these sequences, which the pipeline downloads, but it does not create or update the database.

Configuration
-------------

It is designed to run on the EBI cluster, which as of 2024, uses SLURM.
In theory, it can also run on LSF or any other cluster manager that nextflow supports, but that is untested.

This requires there be a `local.config` file which contains local configuration.
An example configuration file is:

```groovy
// This defaults to '15' but it should be changed to reflect the expected version.
params.version = '15'

// The cluster uses singularity, but the pipeline should work under docker as
// well.
singularity {
  enabled = true
  cacheDir = "$baseDir/singularity"
}

process {
  // This pipeline is containerized and should be run in this container.
  container = 'docker://rfam/genome_download:latest'

  // Currently the SLURM cluster requires setting the expected time and memory
  // cutoffs. These are low ones which work as a default that needs to be
  // raised in some cases.
  executor = 'slurm'
  maxForks = 2000
  time = '1d'
  memory = '1 GB'
}

// These paths are used in generating the new rfam.conf file
params {
  paths {
    production = ''   // This must be filled out
    software = ''     // This must be filled out
  }
}

// Some values which appear to be useful but these could become outdated
executor {
  $slurm {
    queueSize = 10000
    submitRateLimit = '1sec'
  }
  $lsf {
    queueSize = 10000
    submitRateLimit = '1sec'
  }
}

env.PYTHONPATH="$baseDir"

// This should be pointed at ENA's FTP location on the filesystem. This is used
// to determine if some ENA files can be fetched from the filesystem instead of
// using the network.
env.ENA_PATH = ""
```

Note that the `params.paths.production` and `params.paths.software` must be set to generate the production paths for Rfam data and software install respectively.
These values are specific to the cluster setup, which has changed over time.

Running
-------

To run do:

```bash
sbatch run.sh
```

This relies on nextflow fetching the docker image and building a singularity image automatically.

This should produce files in `genomes/` which are the final outputs.
The pipeline will not fail if some genomes are not downloaded.
It is up to the user to validate that the missed genomes are not an issue.

Outputs
-------

This will put all files it downloads under `genomes/` in the base directory and this is not configurable.
The files that are produced are:

- `genomes/failures` - A JSONL formatted file for each job containing all genomes that fail to be downloaded.
    Jobs which have no failed downloads will produce an empty file here.
    To find all failures, use:

    ```sh
    find genomes/failures -type f -size +0 | xargs cat
    ```

- `genomes/fasta` - A directory which contains one directory per job with all downloaded genomes.
    The files are named by the proteome or metagenome that was downloaded

- `genomes/uniprot` - Contains the fetched reference proteome file from uniprot.

- `genomes/metadata` - A series of JSON files of the metadata which can be used to produce the SQL files to update Rfam's mysql database.

- `genomes/rfamseq` - The final files which form rfamseq.

How to use the outputs
----------------------

- Validate that genomes which were not downloaded are acceptable both in quantity and specifically which ones.
- Assuming nothing important, or too many things are missed, create a new rfamseq directory and copy the files in `genomes/rfamseq` into there. This directory must be owned by the Rfam production user.
- Test scanning the new database with a small (in terms of secondary structure size) family.
- Update the mysql database with the new genome metadata.

What this does
--------------

The general task is to fetch the list of [reference proteomes](https://www.uniprot.org/help/reference_proteome) from UniProt.
Generally these genomes are GCA or GCF accessions and can be fetched easily from NCBI/ENA.
The pipeline prefers to fetch from NCBI whenever possible as they have most GCA/GCF sequences organized in a simple FTP structure.

During downloading this will always try to fetch the active version of a genome if possible.
It has been observed that sometimes the specific version of a genome, eg `GCA_N.1` has been suppressed or replaced by a newer version, eg `GCA_N.2`.
This will try to upgrade to the newer ones and log the failure if not possible.
In a few rare cases the genome is suppressed and there is no newer version, in this case nothing is downloaded and the failure is logged.

Unlike previous version of the pipeline, this does not deal with [components](https://www.uniprot.org/help/proteome_component).
The reason is deal with these can be very complicated as it may involve resolving [Whole Genome Sets](https://www.ncbi.nlm.nih.gov/genbank/wgs/) (WGS), which is error-prone.
By removing this task and simply using the GCA/GCF accession the pipeline runs faster, is more reliable and simpler to understand.

Development notes
-----------------

This is a python program which is managed by [`poetry`](https://python-poetry.org/).
An attempt has been made to have the `./pyproject.toml` respect current python standards, but this has not been validated.
If using a package manager that respects the standards, please validate the file first.

All code should be formatted with `black` and `isort`.

There are few external dependencies, which are documented in the Dockerfile.
The docker image is not built automatically, on OSX it can be done with:

```sh
$ docker buildx build -t rfam/genome_download --platform linux/amd64 .
$ docker push
```
