#!/usr/bin/env python

"""
This script fills in any gaps (i.e. no call for whole chromosomes) in the Sclust CNV profile with
a non-value placeholder for downstream analysis
"""
import sys
import re

class SclustProfile(object):

    def __init__(self, profile_line):
        """Parse Sclust profile file to construct object"""
        chrom, start, end, profile_field = profile_line.rstrip().split('\t')
        self.chrom = chrom
        self.start  = start
        self.end = end
        self.profile_field = profile_field
        self.segment_id = self.chrom + ":" + self.start + "-" + self.end

input_args = sys.argv

full_profile_dict = {}
profile_id_to_chrom_dict = {}
with open(input_args[1]) as sclust_file:
    for line in sclust_file:
        sclust_profile_obj = SclustProfile(line)

        full_profile_dict[sclust_profile_obj.segment_id] = (sclust_profile_obj.chrom, sclust_profile_obj.start,
                                                            sclust_profile_obj.end, sclust_profile_obj.profile_field)
        profile_id_to_chrom_dict[sclust_profile_obj.segment_id] = sclust_profile_obj.chrom

chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
               "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16",
               "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]

# Chromosome end coordinates derived from the reference FASTA index file
with open(input_args[2]) as ref_genome_index:
    ref_end_coordinates = ref_genome_index.readlines()

# Creat dictonary of chromosome names and placeholder output
placeholder_dict = {}
for chrom in chromosomes:
    placeholder_dict[chrom] = [chrom, "1", ref_end_coordinates[chromosomes.index(chrom)].rstrip(), "."]

# Loop through all chromosomes and check if present in Sclust profile, if not add the placeholder
for chrom in chromosomes:
    if chrom in profile_id_to_chrom_dict.values():

        # Grab the full profile key for each match and print
        for id_key, chrom_value in profile_id_to_chrom_dict.items():
            if chrom_value == chrom:
                print('\t'.join(full_profile_dict[id_key]))
    else:
        # Print the placeholder if chrom is missing
        print('\t'.join(placeholder_dict[chrom]))
