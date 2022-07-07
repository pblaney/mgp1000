#!/usr/bin/env bash
# Download all the necessary files for use in Battenberg

# Function
battenbergRefDownload() {

	# GC Correction
	wget -O GC_correction_hg38_chr.zip --retry 10 "https://www.dropbox.com/sh/bize1n830t0mgzb/AADQD4DTJOF75YmhBDDoQ9nla/GC_correction_hg38?dl=0&lst=" && \
	unzip GC_correction_hg38_chr.zip && \
	mv 1000G_GC_chr*.txt.gz GC_correction_hg38/ && \
	rm GC_correction_hg38_chr.zip
}

#Function call
battenbergRefDownload