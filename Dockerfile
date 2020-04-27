######################################################
#     MGP1000 Bioinformatics Pipeline Dockerfile     #
#                                                    #
#		Build Nextflow WGS analysis pipeline         #
#			container using Ubuntu 16.04             #
######################################################

# Set Linux OS host image 
FROM ubuntu:16.04

# File Author/ Maintainer
MAINTAINER "Patrick Blaney <patrick.blaney@nyulangone.org>"

# Basic OS set up
RUN apt-get update && apt-get install -y \
	git \
	build-essential

# Pull down current version of MGP1000 repository from GitHub
RUN git clone --depth 1 https://github.com/pblaney/mgp1000.git

# Set working directory
WORKDIR mgp1000/

# Pipeline setup using Makefile
RUN make install-nextflow

# test
RUN nextflow -version
