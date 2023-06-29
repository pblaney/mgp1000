#!/usr/bin/env bash
# This script parses Battenberg output CNV profile to extract complete clonal
# total copy number and major/minor alleles per segment

inputBattenbergBed=$1

finalCnvBed=$3

finalAllelesBed=$4

# Parse the CNV file, first remove header for simplest confroming BED format
grep -v 'seqnames' "${inputBattenbergBed}" \
| \
cut -f 1-3,9 > "${finalCnvBed}"

grep -v 'seqnames' "${inputBattenbergBed}" \
| \
awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$10"/"$11}' > "${finalAllelesBed}"