# Myeloma Genome Project 1000
Comprehensive bioinformatics pipeline for the large-scale collaborative analysis of Multiple Myeloma genomes in an effort to deliniate the broad spectrum of somatic events that occur within them.

## Pipeline Overview
In order to analyze hundreds of matched tumor/normal whole-genome samples across multiple data centers in a consistent manner, a pipeline was created that leverages the workflow management, portability, and reproducibility of Nextflow in conjuction with Docker and Singularity.

<img src="https://github.com/pblaney/mgp1000/blob/master/MGP1000Pipeline.png" width="900">

## Running the Pipeline
The pipeline was developed to be run on various HPCs without concern of environment incompatabilities, version issues, or missing dependencies.
