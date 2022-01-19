#!/usr/bin/env bash
# This script parses a filtered DELLY VCF to split the interchromosomal breakend records to 
# create mate record in VCF format

delly_processed_vcf=$1


for id in `grep 'SVTYPE=BND' $delly_processed_vcf | cut -f 3`
    do
        mate_chrom=$(grep $id $delly_processed_vcf | grep -oE 'CHR2=chr.*?;' | sed 's|CHR2=||' | sed 's|;||')

        mate_pos=$(grep $id $delly_processed_vcf | grep -oE 'POS2=\d*?;' | sed 's|POS2=||' | sed 's|;||')

        mate_id="${id}_mate"
        mate_ref="N"
        mate_alt="<BND>"
        mate_qual="."
        mate_filter="PASS"

        mate_info_svtype="SVTYPE=BND;"

        chr2=$(grep $id $delly_processed_vcf | cut -f 1)
        mate_info_chr2="CHR2=${chr2};"

        pos2=$(grep $id $delly_processed_vcf | cut -f 2)
        mate_info_pos2="POS2=${pos2}"

        mate_format=$(grep $id $delly_processed_vcf | cut -f 9)

        mate_sample=$(grep $id $delly_processed_vcf | cut -f 10)

        echo -e "${mate_chrom}\t${mate_pos}\t${mate_id}\t${mate_ref}\t${mate_alt}\t${mate_qual}\t${mate_filter}\t${mate_info_svtype}${mate_info_chr2}${mate_info_pos2}\t${mate_format}\t${mate_sample}"
    done

