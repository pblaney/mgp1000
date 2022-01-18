#!/usr/bin/env bash
# This script parses a filtered SvABA VCF to find the mate records of interchromosomal breakends then
# appends the input VCF with the cooresponding missing record

svaba_processed_vcf=$1

svaba_unprocessed_vcf=$2

missing_mate_records=$3

for id in `grep -v '^#' $svaba_processed_vcf | cut -f 3 | sed 's|:.||' | sort -k1,1n | uniq`
	do
		have_mate=$(grep -c $id $svaba_processed_vcf)

		if [[ $have_mate == 1 ]]; then
			
			chrom1=$(grep $id $svaba_processed_vcf | cut -f 1)

			chrom2=$(grep $id $svaba_processed_vcf | cut -f 5 | grep -oE 'chr.*?:' | sed 's|:||')

			if [[ $chrom1 != $chrom2 ]]; then
				
				sub_id=$(grep $id $svaba_processed_vcf | cut -f 3 | sed 's|.*:||')

				if [[ $sub_id == 1 ]]; then

					zgrep "${id}:2\t" $svaba_unprocessed_vcf >> $missing_mate_records
				else
					
					zgrep "${id}:1\t" $svaba_unprocessed_vcf >> $missing_mate_records
				fi

			fi

		fi
	done
