SHELL:=/bin/bash

# Set default make call to do nothing
none:

###############################################################################

# Install Git LFS
install-gitlfs-linuxamd64:
	wget https://github.com/git-lfs/git-lfs/releases/download/v3.0.2/git-lfs-linux-amd64-v3.0.2.tar.gz
	tar -zxvf git-lfs-linux-amd64-v3.0.2.tar.gz --exclude=README.md --exclude=CHANGELOG.md --exclude=install.sh
	rm -rf man/
	git lfs install
	rm git-lfs-linux-amd64-v3.0.2.tar.gz

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

# Remove logs/pid/reports/trace/slurmsub files
quick-clean:
	rm -f .nextflow.log*
	rm -f .nextflow.pid*
	rm -f nextflow_report.*.html*
	rm -f timeline_report.*.html*
	rm -f trace.*.txt*
	rm -f slurmsub.*.*.*

# Clean up the pipeline directory after a successful run to prep for new run
clean-postrun: quick-clean
	rm -rf work/*

# Completely scrub the pipeline directory
clean-all: clean-postrun
	rm -rf output/*

###############################################################################
