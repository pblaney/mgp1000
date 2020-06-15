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
