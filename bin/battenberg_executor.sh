#!/usr/bin/env bash
# Execute the Battenberg process

# Grab input paramters as the variables for Battenberg command options
tumorId=$1

normalId=$2

tumorBam=$3

normalBam=$4

sex=$5

outputDir=$6

cpus=$7

minDepth=$8

# Update the path to the reference files for use
sed -E 's|\$\{REF_PATH\}|'${PWD}'/battenberg_reference|g' battenberg_reference/impute_info.txt > battenberg_reference/imputation/impute_info.txt

# Second, update the path of the Beagle5 base directory
cat /opt/battenberg/inst/example/battenberg_wgs.R | \
sed 's|BEAGLE_BASEDIR = \".*|BEAGLE_BASEDIR = \"'${PWD}'/battenberg_reference\"|' | \
sed 's|1kg.phase3.v5a_GRCh38nounref_allele_index_|1kg.phase3.v5a_GRCh38nounref_allele_index_chr|' | \
sed 's|1kg.phase3.v5a_GRCh38nounref_loci_|1kg.phase3.v5a_GRCh38nounref_loci_chr|' | \
sed 's|1000G_GC_|1000G_GC_chr|' | \
sed 's|1000G_RT_|1000G_RT_chr|' | \
sed 's|CHROM_COORD_FILE = \".*|CHROM_COORD_FILE = \"/opt/battenberg/chromosome_coordinates_hg38.txt\"|' | \
sed 's|MIN_NORMAL_DEPTH = 10|MIN_NORMAL_DEPTH = '${minDepth}'|' > battenberg_wgs.R

# Make user-defined directory for output
mkdir -p "${outputDir}"

# Execute Battenberg command
Rscript battenberg_wgs.R -a "paired" -t "${tumorId}" -n "${normalId}" --tb "${PWD}/${tumorBam}" --nb "${PWD}/${normalBam}" --sex "${sex}" -o "${outputDir}" --cpu ${cpus} -g "hg38"
