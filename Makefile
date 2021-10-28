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
	gunzip -q references/hg38/Homo_sapiens_assembly38.fasta.gz
	gunzip -q references/hg38/Homo_sapiens_assembly38.fasta.64.bwt.gz
	gunzip -q references/hg38/Homo_sapiens_assembly38.fasta.64.sa.gz
	gunzip -q references/hg38/Homo_sapiens_assembly38_autosome_sex_chroms/*.fa.gz
	mkdir -p input
	mkdir -p input/preprocessedBams
	mkdir -p logs

###############################################################################

# Save the necessary output files and clean the directories of any unneeded files after 
# successful completion of the steps of the pipeline
preprocessing-completion:
	mkdir -p logs/preprocessing
	mv nextflow_report.*.html logs/preprocessing
	mv timeline_report.*.html logs/preprocessing
	mv trace.*.txt logs/preprocessing

germline-completion:
	mkdir -p logs/germline
	mv nextflow_report.*.html logs/germline
	mv timeline_report.*.html logs/germline
	mv trace.*.txt logs/germline

somatic-completion:
	mkdir -p logs/somatic
	mv nextflow_report.*.html logs/somatic
	mv timeline_report.*.html logs/somatic
	mv trace.*.txt logs/somatic

###############################################################################

# Remove logs/pid/reports/trace files
quick-clean:
	rm -f .nextflow.log*
	rm -f .nextflow.pid*
	rm -f timeline_report.*.html*
	rm -f nextflow_report.*.html*
	rm -f trace.*.txt*

# Clean up the pipeline directory after a successful run to prep for new run
clean-postrun: quick-clean
	rm -rf work/*

# Completely scrub the pipeline directory
clean-all:
	rm -rf work
	rm -rf output
	rm -f .nextflow.log*
	rm -f .nextflow.pid*
	rm -f timeline_report.*.html*
	rm -f nextflow_report.*.html*
	rm -f trace.*.txt*

###############################################################################
