#!/usr/bin/env bash
# This script parses FACETS output CNV profile to extract complete clonal
# total copy number and major/minor alleles per segment

inputFacetsProfile=$1

sampleId=$2

finalCnvBed=$3

finalAllelesBed=$4

# Parse the CNV file
grep -v 'chrom' ${inputFacetsProfile} \
| \
awk 'BEGIN {OFS="\t"} {print "chr"$1,$10,$11,$13,$14}' \
| \
sed 's|^chr23|chrX|' > "${sampleId}.facets.cnv.simplified.bed"

# Easy grab of final CNV bed
cut -f 1-4 "${sampleId}.facets.cnv.simplified.bed" > "${finalCnvBed}"

# Read each segment and output the infered major/minor alleles of the NA calls
# The inference will assume heterozygosity with a single minor allele
while read -r cnvSegment
	do
    	# Check if the minor allele column contains an NA instead of a value
    	minorAllele=$(echo ${cnvSegment} | cut -d ' ' -f 5)
    	totalCn=$(echo ${cnvSegment} | cut -d ' ' -f 4)
    	hetAssumptionCn=$(( $totalCn - 1 ))

    	if [[ ${minorAllele} == "NA" ]]; then

    		awk -v inferredCn="$hetAssumptionCn" 'BEGIN {OFS="\t"} {print $1,$2,$3,inferredCn"/1"}' <(echo ${cnvSegment})

    	else
    		# Some segments report the minor allele as more than the major (total_cn - minor_allele)
    		majorAllele=$(( $totalCn - $minorAllele ))

    		if [[ $(echo $minorAllele) > $(echo $majorAllele) ]]; then
    			awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$5"/"$4-$5}' <(echo ${cnvSegment})
    		else
    			awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4-$5"/"$5}' <(echo ${cnvSegment})
    		fi

    	fi
    done < "${sampleId}.facets.cnv.simplified.bed" 1> "${finalAllelesBed}"
