SHELL:=/bin/bash

# Set default make call to do nothing
none:

###############################################################################

# Install/update Nextflow 
./nextflow:
	curl -fsSL get.nextflow.io | bash

install-nextflow: ./nextflow

update-nextflow: ./nextflow
	./nextflow self-update

###############################################################################

# Prepare reference genome files
prep-ref:
	gunzip references/hg38/bwa/genome.fa.gz
	gunzip references/hg38/bwa/genome.fa.bwt.gz
	gunzip references/hg38/bwa/genome.fa.sa.gz

# Test pipeline locally with Docker and BAM input files
local-test-run-bam:
	nextflow run main.nf -resume --input_format bam -profile docker 

# Test pipeline locally with Docker and FASTQ input files
local-test-run-fastq:
	nextflow run main.nf -resume --input_format fastq -profile docker 

# Remove all output files
clean-all:
	rm -f .nextflow.log*
	rm -rf work
	rm -rf output
	rm -f timeline_report.html*
	rm -f nextflow_report.html*
	rm -f trace.txt*
