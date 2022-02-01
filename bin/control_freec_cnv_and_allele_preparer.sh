#!/usr/bin/env bash
# Create a complete total copy number and allele genotype BED file from Control-FREEC bedGraph
# output using BEDtools

tumor_normal_sample_id=$1

control_freec_cnv_txt=$2

control_freec_bedgraph=$3

reference_fasta_index=$4

# Prep the Control-FREEC total copy number and alleles from the *controlfreec.cnv.txt file
# Total copy number
grep -v 'chr' "${control_freec_cnv_txt}" \
| \
awk 'BEGIN {OFS="\t"} {print "chr"$1,$2,$3,$4}' \
| \
sort -k1,1V -k2,2n > "${tumor_normal_sample_id}.controlfreec.cnv.bed"

# Alleles
grep -v 'chr' "${control_freec_cnv_txt}" \
| \
awk 'BEGIN {OFS="\t"} {print "chr"$1,$2,$3,$6}' \
| \
sort -k1,1V -k2,2n > "${tumor_normal_sample_id}.controlfreec.alleles.bed"

# To split the Control-FREEC bedGraph into files based on track (CN LOH, gains, losses, CN neutal)

csplit -s -z -f "${tumor_normal_sample_id}".controlfreec.cn.track "${control_freec_bedgraph}" /track/ "{$(($(grep -c 'track' "${control_freec_bedgraph}") - 1))}"

# Loop through track files to determine which contains which specific track information
for cn_track_file in `ls -1 *.controlfreec.cn.track*`;
    do
        if grep -q 'copy neutral LOH' $cn_track_file; then
            grep -v 'track' $cn_track_file | sort -k1,1V -k2,2n > "${tumor_normal_sample_id}".controlfreec.nloh.track.bed
        fi

        if grep -q 'gains' $cn_track_file; then
           grep -v 'track' $cn_track_file | sort -k1,1V -k2,2n > "${tumor_normal_sample_id}".controlfreec.gain.track.bed
        fi

        if grep -q 'losses' $cn_track_file; then
            grep -v 'track' $cn_track_file | sort -k1,1V -k2,2n > "${tumor_normal_sample_id}".controlfreec.loss.track.bed
        fi

        if grep -q 'neutral copy number' $cn_track_file; then
           grep -v 'track' $cn_track_file | sort -k1,1V -k2,2n > "${tumor_normal_sample_id}".controlfreec.neutral.track.bed
        fi
     done

# For converting the the copy number neutral track to a total copy number and allele BED
# First, calculate the average of the copy number column and then round this to a whole value
neutral_cn_value=$(awk '{ total += $4; count++ } END { print total/count }' "${tumor_normal_sample_id}".controlfreec.neutral.track.bed)
rounded_neutral_cn_value=$(printf "%.0f" $neutral_cn_value)

awk -v segment_cn=$rounded_neutral_cn_value 'BEGIN {OFS="\t"} {print $1,$2,$3,segment_cn}' "${tumor_normal_sample_id}".controlfreec.neutral.track.bed > "${tumor_normal_sample_id}".controlfreec.neutral.cnv.mapped.bed

# Then use the rounded total copy number value to determine the possible alleles 
if [[ $rounded_neutral_cn_value == 2 ]]; then
     awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"AB"}' "${tumor_normal_sample_id}".controlfreec.neutral.track.bed > "${tumor_normal_sample_id}".controlfreec.neutral.alleles.mapped.bed
elif [[ $rounded_neutral_cn_value == 3 ]]; then
     awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"AAB"}' "${tumor_normal_sample_id}".controlfreec.neutral.track.bed > "${tumor_normal_sample_id}".controlfreec.neutral.alleles.mapped.bed
elif [[ $rounded_neutral_cn_value == 4 ]]; then
     awk 'BEGIN {OFS="\t"} {print $1,$2,$3,"AABB"}' "${tumor_normal_sample_id}".controlfreec.neutral.track.bed > "${tumor_normal_sample_id}".controlfreec.neutral.alleles.mapped.bed
fi

# For converting the copy number neutral LOH track to a total copy number and allele BED
# Map the copy number values and allele states from BED files derived from the *controlfreec.cnv.txt file
bedtools map \
-a "${tumor_normal_sample_id}".controlfreec.nloh.track.bed \
-b "${tumor_normal_sample_id}.controlfreec.cnv.bed" \
-c 4 \
-o distinct \
-g "${reference_fasta_index}" \
| \
cut -f 1-3,5 \
| \
sort -k1,1V -k2,2n > "${tumor_normal_sample_id}".controlfreec.nloh.cnv.mapped.bed

bedtools map \
-a "${tumor_normal_sample_id}".controlfreec.nloh.track.bed \
-b "${tumor_normal_sample_id}.controlfreec.alleles.bed" \
-c 4 \
-o distinct \
-g "${reference_fasta_index}" \
| \
cut -f 1-3,5 \
| \
sort -k1,1V -k2,2n > "${tumor_normal_sample_id}".controlfreec.nloh.alleles.mapped.bed

# For converting the copy number loss track to a total copy numver and allele BED
# Map the copy number values and allele states from BED files derived from the *controlfreec.cnv.txt file
# However, to account for situations with a '.' as the total copy number, extract these lines separately 
# and infer the copy number by rounding the bedGraph value
bedtools map \
-a "${tumor_normal_sample_id}".controlfreec.loss.track.bed \
-b "${tumor_normal_sample_id}.controlfreec.cnv.bed" \
-c 4 \
-o distinct \
-g "${reference_fasta_index}" \
| \
grep -vP "\.$" \
| \
cut -f 1-3,5 > "${tumor_normal_sample_id}".controlfreec.loss.cnv.mapped.intermediate1.bed

bedtools map \
-a "${tumor_normal_sample_id}".controlfreec.loss.track.bed \
-b "${tumor_normal_sample_id}.controlfreec.cnv.bed" \
-c 4 \
-o distinct \
-g "${reference_fasta_index}" \
| \
grep -P "\.$" \
| \
awk '{printf "%s\t%s\t%s\t%.0f\n", $1,$2,$3,$4}' > "${tumor_normal_sample_id}".controlfreec.loss.cnv.mapped.intermediate2.bed

# Merge the two intermediate loss total copy number BED files
cat "${tumor_normal_sample_id}".controlfreec.loss.cnv.mapped.intermediate1.bed "${tumor_normal_sample_id}".controlfreec.loss.cnv.mapped.intermediate2.bed \
| \
sort -k1,1V -k2,2n > "${tumor_normal_sample_id}".controlfreec.loss.cnv.mapped.bed

# Also, to account for '.' as alleles for segments with single copy loss, infer the alleles using
# the mapped total copy number
bedtools map \
-a "${tumor_normal_sample_id}".controlfreec.loss.cnv.mapped.bed \
-b "${tumor_normal_sample_id}.controlfreec.alleles.bed" \
-c 4 \
-o distinct \
-g "${reference_fasta_index}" \
| \
sed 's|1\t\.|1\tA|' \
| \
cut -f 1-3,5 \
| \
sort -k1,1V -k2,2n > "${tumor_normal_sample_id}".controlfreec.loss.alleles.mapped.bed

# For converting the copy number gain track to a total copy numver and allele BED
# Map the copy number values and allele states from BED files derived from the *controlfreec.cnv.txt file
# However, to account for situations with a '.' as the total copy number, extract these lines separately 
# and infer the copy number by rounding the bedGraph value
bedtools map \
-a "${tumor_normal_sample_id}".controlfreec.gain.track.bed \
-b "${tumor_normal_sample_id}.controlfreec.cnv.bed" \
-c 4 \
-o distinct \
-g "${reference_fasta_index}" \
| \
grep -vP "\.$" \
| \
cut -f 1-3,5 > "${tumor_normal_sample_id}".controlfreec.gain.cnv.mapped.intermediate1.bed

bedtools map \
-a "${tumor_normal_sample_id}".controlfreec.gain.track.bed \
-b "${tumor_normal_sample_id}.controlfreec.cnv.bed" \
-c 4 \
-o distinct \
-g "${reference_fasta_index}" \
| \
grep -P "\.$" \
| \
awk '{printf "%s\t%s\t%s\t%.0f\n", $1,$2,$3,$4}' > "${tumor_normal_sample_id}".controlfreec.gain.cnv.mapped.intermediate2.bed

# Merge the two intermediate gain total copy number BED files
cat "${tumor_normal_sample_id}".controlfreec.gain.cnv.mapped.intermediate1.bed "${tumor_normal_sample_id}".controlfreec.gain.cnv.mapped.intermediate2.bed \
| \
sort -k1,1V -k2,2n > "${tumor_normal_sample_id}".controlfreec.gain.cnv.mapped.bed

bedtools map \
-a "${tumor_normal_sample_id}".controlfreec.gain.track.bed \
-b "${tumor_normal_sample_id}.controlfreec.alleles.bed" \
-c 4 \
-o distinct \
-g "${reference_fasta_index}" \
| \
cut -f 1-3,5 \
| \
sort -k1,1V -k2,2n > "${tumor_normal_sample_id}".controlfreec.gain.alleles.mapped.bed

# Now that all total copy number and allele BEDs are in same format, merge them for mapping
cat "${tumor_normal_sample_id}".controlfreec.neutral.cnv.mapped.bed "${tumor_normal_sample_id}".controlfreec.nloh.cnv.mapped.bed "${tumor_normal_sample_id}".controlfreec.loss.cnv.mapped.bed "${tumor_normal_sample_id}".controlfreec.gain.cnv.mapped.bed \
| \
sort -k1,1V -k2,2n > "${tumor_normal_sample_id}".controlfreec.complete.cnv.mapped.bed

cat "${tumor_normal_sample_id}".controlfreec.neutral.alleles.mapped.bed "${tumor_normal_sample_id}".controlfreec.nloh.alleles.mapped.bed "${tumor_normal_sample_id}".controlfreec.loss.alleles.mapped.bed "${tumor_normal_sample_id}".controlfreec.gain.alleles.mapped.bed \
| \
sort -k1,1V -k2,2n > "${tumor_normal_sample_id}".controlfreec.complete.alleles.mapped.bed

# Finally, map the complete total copy number and allele BEDs together for final adjustment within python script
bedtools map \
-a "${tumor_normal_sample_id}".controlfreec.complete.cnv.mapped.bed \
-b "${tumor_normal_sample_id}".controlfreec.complete.alleles.mapped.bed \
-c 4 \
-o distinct > "${tumor_normal_sample_id}".controlfreec.complete.cnv.alleles.merged.bed
