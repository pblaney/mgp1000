#!/usr/bin/env bash
# Execute the fragCounter process

# Grab input parameters as the variables for fragCounter command options
inputBam=$1

gcMappabilityDir=$2

isPaired=$3

subChrString=$4

outputDir=$5

# Remove she-bang header to allow for execution as vanilla .R script
grep -v 'Rscript' /opt/toolshed/R/library/fragCounter/extdata/frag > fragCounter.R
chmod 755 fragCounter.R

# Make user-defined directory for output
mkdir -p "${outputDir}"

# Execute fragCounter command
Rscript --vanilla fragCounter.R \
--bam "${inputBam}" \
--window 200 \
--gcmapdir "${gcMappabilityDir}" \
--paired "${isPaired}" \
--chrsub "${subChrString}" \
--outdir "${outputDir}"
