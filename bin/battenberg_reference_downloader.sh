#!/usr/bin/env bash
# Download all the necessary files for use in Battenberg

# Debugging settings
set -euo pipefail

echo "###########################################################"
echo "#          Battenberg Reference Data Downloader           #"
echo "###########################################################"
echo 

# Pull whole zip file
echo "Beginning Download ..."
wget -O battenberg_reference_hg38.zip "https://ora.ox.ac.uk/objects/uuid:08e24957-7e76-438a-bd38-66c48008cf52/files/rt435gd52w" && \
7zz x battenberg_reference_hg38.zip
rm battenberg_reference_hg38.zip
echo "Done"

# GC Correction
echo "Prepping GC Correction ..."
unzip GC_correction_hg38.zip -x
rm GC_correction_hg38.zip
rm GC_correction_hg38_chr.zip
rm GC_correction_hg38/1000G_GC_chr23.txt.gz
echo "Done"

# RT Correction
echo "Prepping RT Correction ..."
unzip RT_correction_hg38.zip -x
rm RT_correction_hg38.zip
rm RT_correction_hg38_chr.zip
rm RT_correction_hg38/1000G_RT_chr23.txt.gz
echo "Done"

# Shapeit2
echo "Prepping Shapeit2 ..."
unzip shapeit2_chr.zip -x
rm shapeit2_chr.zip
rm shapeit2.zip
rm shapeit2/ALL.v1a.shapeit2_integrated_chr23*
echo "Done"

# Impute
echo "Prepping Imputation and Impute info ..."
unzip imputation_chr.zip -x
rm imputation_chr.zip
rm imputation.zip
echo "Done"

# 1000G
echo "Prepping 1000 Genome Loci ..."
unzip 1000G_loci_hg38_chr.zip -x
rm 1000G_loci_hg38_chr.zip
rm 1000G_loci_hg38.zip
rm 1000G_loci_hg38/1kg.phase3.v5a_GRCh38nounref_*chr23*
echo "Done"

# Probloci
echo "Prepping Problem Loci ..."
unzip probloci.zip -x
rm probloci_chr.zip
rm probloci.zip
echo "Done"

# Beagle5
echo "Prepping Beagle5 ..."
unzip beagle_chr.zip -x
rm beagle.zip
rm beagle_chr.zip
rm beagle/chr23.1kg.phase3.v5a_GRCh38nounref.vcf.gz
rm beagle/plink.chr23*
rm beagle/plink.chr23*
# Need to add 'chr' notation to VCFs
mkdir -p temp/
mv beagle/chr*.1kg.phase3.v5a_GRCh38nounref.vcf.gz temp/
for i in {1..22} X; do zcat "temp/chr${i}.1kg.phase3.v5a_GRCh38nounref.vcf.gz" | sed -E 's|^'${i}'|chr'${i}'|' | gzip > "beagle/chr${i}.1kg.phase3.v5a_GRCh38nounref.vcf.gz"; done
rm -rf temp/
mv beagle/ beagle5/
echo "Done"

# Final clean up
rm sha1sums.txt
rm chromosome_coordinates_hg38_chr.txt
rm chromosome_coordinates_hg38.txt

echo 
echo "Completed Battenberg reference data preparation"
echo 
