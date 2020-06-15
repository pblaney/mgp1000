# Myeloma Genome Project 1000
Comprehensive bioinformatics pipeline for the large-scale collaborative analysis of Multiple Myeloma genomes in an effort to deliniate the broad spectrum of somatic events

## Pipeline Overview
In order to analyze over one thousand matched tumor/normal whole-genome samples across multiple data centers in a consistent manner, a pipeline was created that leverages the workflow management, portability, and reproducibility of [Nextflow](http://www.nextflow.io/) in conjuction with [Singularity](https://sylabs.io/docs/).

The entire pipeline is divided into 3 steps: Preprocessing, Germline Variant Analysis, and Somatic Variant Analysis
This compartmentalizes the workflow and provides significant completion checkpoints which is effective for large-scale batch processing. 

<img src="https://github.com/pblaney/mgp1000/blob/master/MGP1000Pipeline.png" width="900">

## Running the Pipeline
The pipeline was developed to be run on various HPCs without concern of environment incompatabilities, version issues, or missing dependencies. However, there are a few assumptions regarding initial setup of the pipeline but the required software should be readily available in nearly all HPC environments.
* Git
* GNU Utilities
* Java 8 (or later)
* Singularity v3.1

# Clone GitHub Repository
The first step in the deployment process is to clone the MGP1000 GitHub repository to a location on your HPC that is large enough to hold the input/output data, like a scratch directory, and has access to the job scheduling software, such as Slurm or SGE.
```
$ cd <scratch dir>

$ git clone https://github.com/pblaney/mgp1000.git
### Example output ###
# Cloning into 'mgp1000'...
# remote: Enumerating objects: 83, done.
# remote: Counting objects: 100% (83/83), done.
# remote: Compressing objects: 100% (67/67), done.
# remote: Total 300 (delta 29), reused 64 (delta 13), pack-reused 217
# Receiving objects: 100% (300/300), 412.90 KiB | 13.76 MiB/s, done.
# Resolving deltas: 100% (155/155), done.

$ cd mgp1000/
```

# Install Nextflow
This series of `make` commands will install Nextflow, test, and, optionally, update the current Nextflow installation. First, check for what current version of Java is available to the current environment.
```
$ java -version
### Example output ###
# openjdk version "1.8.0_131"
# OpenJDK Runtime Environment (build 1.8.0_131-b12)
# OpenJDK 64-Bit Server VM (build 25.131-b12, mixed mode)

$ make install-nextflow
### Example output ###
# curl -fsSL get.nextflow.io | bash
# CAPSULE: Downloading dependency .....
# ....
# ....
# 	  N E X T F L O W
#     version 20.04.1 build 5335
#     created 03-05-2020 19:37 UTC (15:37 EDT)
#     cite doi:10.1038/nbt.3820
#     http://nextflow.io

$ make nextflow-test
### Example output ###

```

# Prepare the Pipeline for Usage
Due to size, certain reference genome files are GNU zipped so these `make` commands will unzip them for use in the pipeline. Additionally, an `input` directory is created for staging all input BAM or FASTQ files.
```
$ make prep-pipeline
```
