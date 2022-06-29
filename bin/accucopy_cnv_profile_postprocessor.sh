#!/usr/bin/env bash
# This script parses an Accucopy output CNV profile to extract complete clonal and subclonal
# total copy number and major/minor alleles per segment

inputAccucopyProfile=$1

sampleId=$2

subclonalOutput=$3

finalCnvBed=$4

finalAllelesBed=$5

# First, separately extract all the clonal and subclonal calls
grep -Ev '^chr|^#' ${inputAccucopyProfile} \
| \
grep -v 'NA' \
| \
awk 'BEGIN {OFS="\t"} {print "chr"$1,$2,$3,$4,$5,$4-$5}' \
| \
sort -k1,1V -k2,2n > "${sampleId}.accucopy.clonal.bed"

grep -Ev '^chr|^#' ${inputAccucopyProfile} \
| \
grep 'NA' \
| \
awk 'BEGIN {OFS="\t"} {print "chr"$1,$2,$3,$4}' \
| \
sort -k1,1V -k2,2n > "${subclonalOutput}"

# Round the subclonal segments to generate the likely clonal total copy number per segment
paste <(cut -f 1-3 "${subclonalOutput}") <(printf "%.0f\n" $(cut -f 4 "${subclonalOutput}")) \
| \
sed 's|\-1|0|' > "${sampleId}.accucopy.subclonal.rounded.bed"

# Read each segment and output the infered major/minor alleles of the subclonal calls
# The inference will assume heterozygosity with a single minor allele
while read -r subclonalSegment
	do
    	roundedCn=$(echo ${subclonalSegment} | cut -d ' ' -f 4)
    	hetAssumptionCn=(( $roundedCn - 1 ))

    	if [[ ${roundedCn} == 0 ]]; then
      		awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,"0","0"}' <(echo ${subclonalSegment})
    	elif [[ ${roundedCn} == 1 ]]; then
      		awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,"1","0"}' <(echo ${subclonalSegment})
    	else
      		awk -v inferredCn="$hetAssumptionCn" 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,inferredCn,"1"}' <(echo ${subclonalSegment})
    	fi
  	done < "${sampleId}.accucopy.subclonal.rounded.bed" 1> "${sampleId}.accucopy.subclonal.rounded.alleles.bed"

# Combine the clonal calls and the rounded subclonal calls with inferred major/minor alleles
cat "${sampleId}.accucopy.clonal.bed" "${sampleId}.accucopy.subclonal.rounded.alleles.bed" \
| \
sort -k1,1V -k2,2n > "${sampleId}.accucopy.cnv.alleles.bed"

# Subset the final CNV profile bed into separate BED files for the total copy number and major/minor alleles
# so they can be incorporated in the consensus mechanism
awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$4}' "${sampleId}.accucopy.cnv.alleles.bed" > "${finalCnvBed}"

awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$5"/"$6}' "${sampleId}.accucopy.cnv.alleles.bed" \
| \
sed 's|\.\/\.|.|' > "${finalAllelesBed}"
