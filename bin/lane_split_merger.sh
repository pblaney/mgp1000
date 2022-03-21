#!/usr/bin/env bash
# Given an input directory, merge all lane-split FASTQ files per R1/R2 sets

inputDirectory=$1

outputDirectory=$2

mkdir -p "${outputDirectory}"
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
	find -L "${inputDirectory}" -maxdepth 8 -type f -name "${P}_00*.fastq.gz" -exec cat '{}' ';' > "${P}.merged.fastq.gz"
	mv "${P}.merged.fastq.gz" "${outputDirectory}"
done
