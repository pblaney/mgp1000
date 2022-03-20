#!/usr/bin/env bash
# Given an input directory, merge all lane-split FASTQ files per R1/R2 sets

inputDirectory=$1

outputDirectory=$2

mkdir -p "${outputDirectory}"
mkdir -p tmp/
ls -1 ${inputDirectory}/*.fastq.gz \
| \
while read F
do 
	basename $F | sed -E 's|_00[1234].fastq.gz||'
done \
| \
sort \
| \
uniq \
| \
while read P
do
	find -L "${inputDirectory}" -maxdepth 8 -type f -name "${P}_00*.fastq.gz" -exec cat '{}' ';' > "${P}.unsorted.merged.fastq.gz"
	gunzip "${P}.unsorted.merged.fastq.gz"
	fastq-sort --temporary-directory tmp/ --idn "${P}.unsorted.merged.fastq" | bgzip > "${P}.merged.fastq.gz"
	mv "${P}.merged.fastq.gz" "${outputDirectory}"
	rm "${P}.unsorted.merged.fastq"
done
