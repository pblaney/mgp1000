SHELL:=/bin/bash

# Set default make call to do nothing
none:

###############################################################################

# Install/update Nextflow 
./nextflow:
	curl -fsSL get.nextflow.io | bash

install-nextflow: ./nextflow

nextflow-test: ./nextflow
	./nextflow run hello

update-nextflow: ./nextflow
	./nextflow self-update

###############################################################################

# Prepare reference genome files and create input directory
prep-pipeline:
	gunzip references/hg38/Homo_sapiens_assembly38.fasta.gz
	gunzip references/hg38/Homo_sapiens_assembly38.fasta.64.bwt.gz
	gunzip references/hg38/Homo_sapiens_assembly38.fasta.64.sa.gz
	mkdir -p input
	mkdir -p input/preprocessedBams
	mkdir -p logs

###############################################################################

# Test run the pipeline within production environment to validate setup
test-run:
	nextflow run preprocessing.nf -bg -resume --input_format fastq -profile preprocessing

###############################################################################

# Save the necessary output files and clean the directory of any unneeded files after 
# successful completion of the Preprocessing step
preprocessing-completion:
	mkdir -p logs/preprocessing
	mv nextflow_report.html logs/preprocessing
	mv timeline_report.html logs/preprocessing
	mv trace.txt logs/preprocessing
	mv output/preprocessing/finalPreprocessedBams/* input/preprocessedBams
	rm -rf work/*

###############################################################################

# Test Preprocessing step locally with Docker and BAM or FASTQ input files
dev-preprocessing-bam:
	nextflow run preprocessing.nf -bg -resume --input_format bam -profile dev_preprocessing

dev-preprocessing-fastq:
	nextflow run preprocessing.nf -bg -resume --input_format fastq -profile dev_preprocessing

dev-preprocessing-bigbam:
	nextflow run preprocessing.nf -bg --input_format bam --skip_to_qc yes -profile dev_preprocessing

# Test Germline step locally with Docker
dev-germline:
	nextflow run germline.nf -bg -resume -profile dev_germline

###############################################################################

# Remove logs/pid/reports/trace files
quick-clean:
	rm -f .nextflow.log*
	rm -f .nextflow.pid*
	rm -f timeline_report.html*
	rm -f nextflow_report.html*
	rm -f trace.txt*

# Completely scrub pipeline output files
clean-all:
	rm -rf work
	rm -rf output
	rm -f .nextflow.log*
	rm -f .nextflow.pid*
	rm -f timeline_report.html*
	rm -f nextflow_report.html*
	rm -f trace.txt*

# Scrub Preprocessing step output files
clean-preprocessing:
	rm -rf work/*
	rm -rf output/preprocessing/*
	rm -f .nextflow.log*
	rm -f .nextflow.pid*
	rm -f timeline_report.html*
	rm -f nextflow_report.html*
	rm -f trace.txt*

# Scrub Germline step output files
clean-germline:
	rm -rf work/*
	rm -rf output/germline/*
	rm -f .nextflow.log*
	rm -f .nextflow.pid*
	rm -f timeline_report.html*
	rm -f nextflow_report.html*
	rm -f trace.txt*

###############################################################################