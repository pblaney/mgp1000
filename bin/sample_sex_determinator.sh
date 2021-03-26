#!/usr/bin/env bash
# Determine the sex of the sample for use in Battenberg
# This code was adopted from the Bioinformatic support for Cancer, Ageing and Somatic Mutation group
# at the Wellcome Trust Sanger Institute's implementation of Battenberg
# https://github.com/cancerit/cgpBattenberg/blob/dev/perl/lib/Sanger/CGP/Battenberg/Implement.pm

allele_counter_output=$1

# Parse the AlleleCounter output file for depth > 5 at any of the gender identification loci
depth_at_loci=$(cut -f 7 "$allele_counter_output" | tail -n 4)

gender_determinator() {
	gender="female XX"
	for depth in $depth_at_loci; do
		if [[ $depth -gt 5 ]]; then
			gender="male XY"
		fi
	done
	echo "$gender"
}

gender_determinator
