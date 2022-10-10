#!/usr/bin/env bash
# Given an input directory, merge all lane-split FASTQ files per R1/R2 sets

inputDirectory=$1

workDirectory=$2

outputDirectory=$3

mkdir -p "${outputDirectory}"
ls -1 ${workDirectory}/*.f*q.gz \
| \
while read fastqPath
do
	# Find base naming convention
	basename $fastqPath | sed -E 's|\.f.*q.gz||'
done \
| \
sort \
| \
while read fastqBasename
do
	# Test different lane split conventions
	if [[ ${fastqBasename} =~ _L00[1-9]_R[12]+ ]]; then

		echo ${fastqBasename} | sed -E 's|_L00[1-9]_R|_L00\*_R|'
	
	elif [[ ${fastqBasename} =~ _00[1-9]_R[12]+ ]]; then

		echo ${fastqBasename} | sed -E 's|_00[1-9]_R|_00\*_R|'

	elif [[ ${fastqBasename} =~ _R[12]_00[1-9] ]]; then

		echo ${fastqBasename} | sed -E 's|_00[1-9]$|_00\*|'
	fi
done \
| \
sort \
| \
uniq \
| \
while read laneSplitRegex
do
	# Concat lane split FASTQs
	splitFiles=$(find -L "${inputDirectory}" -maxdepth 8 -type f -name "${laneSplitRegex}.f*q.gz" | sort)

	if [[ ${laneSplitRegex} =~ _L00\*_R ]]; then

		laneMergedFastqBasename=$(echo "$laneSplitRegex" | sed -E 's|_L00\*||')

	elif [[ ${laneSplitRegex} =~ _00\*_R ]]; then

		laneMergedFastqBasename=$(echo "$laneSplitRegex" | sed -E 's|_00\*||')

	elif [[ ${laneSplitRegex} =~ _R[12]_00\* ]]; then

		laneMergedFastqBasename=$(echo "$laneSplitRegex" | sed -E 's|_00\*||')
	fi

	mergeCmd="cat ${splitFiles} > ${laneMergedFastqBasename}.merged.fastq.gz"
	eval ${mergeCmd}
	mv "${laneMergedFastqBasename}.merged.fastq.gz" "${outputDirectory}"
done
