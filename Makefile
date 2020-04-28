SHELL:=/bin/bash

# Set default make call to do nothing
none:

# Install/update Nextflow 
install-nextflow:
	curl -fsSL get.nextflow.io | bash

update-nextflow:
	./nextflow self-update