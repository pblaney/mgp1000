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

# Get name of Beagle5 JAR file
BEAGLEJAR=$(ls "${PWD}/battenberg_reference/beagle5/beagle*.jar" | sed -E 's|\/.*\/||')

# Second, update the path of the Beagle5 base directory
cat /opt/battenberg/inst/example/battenberg_wgs.R | \
sed 's|BEAGLE_BASEDIR = \".*|BEAGLE_BASEDIR = \"'${PWD}'/battenberg_reference\"|' | \
sed 's|beagle.*.jar|'${BEAGLEJAR}'|' | \
sed 's|CHROM_COORD_FILE = \".*|CHROM_COORD_FILE = \"/opt/battenberg/chromosome_coordinates_hg38.txt\"|' > battenberg_wgs.R

# Make user-defined directory for output
mkdir -p "${outputDir}"

# Execute Battenberg command
cmd="
Rscript battenberg_wgs.R \
--analysis_type paried \
--tumourname ${tumorId} \
--normalname ${normalId} \
--tb ${PWD}/${tumorBam} \
--nb ${PWD}/${normalBam} \
--sex ${sex} \
--output ${outputDir} \
--cpu ${cpus} \
--ref_genome_build hg38 \
--min_depth ${minDepth}
"
echo "${cmd}"
eval ${cmd}
