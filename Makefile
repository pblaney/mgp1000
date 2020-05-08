SHELL:=/bin/bash

# Set default make call to do nothing
none:

# Install/update Nextflow 
./nextflow:
	curl -fsSL get.nextflow.io | bash

install-nextflow: ./nextflow

update-nextflow: ./nextflow
	./nextflow self-update

# Test pipeline locally with Docker
local-test-run:
	nextflow run main.nf -profile docker

# Remove all output files
clean-all:
	rm -f .nextflow.log*
	rm -rf work
	rm -rf output
	rm -f timeline_report.html*
	rm -f nextflow_report.html*
	rm -f trace.txt*