#!/usr/bin/env python

"""
This script calculates a forward-reverse score for each site within a VCF in preparation for use in Sclust.
The forward-reverse score is defined as 0 if all reads of the mutation are facing in one direction and
1 if all forward and reverse scores are equally present.
"""
import fileinput

for line in fileinput.input():
    tumor_forward_reverse_read_attributes = line.split('\t')
    chromosome = tumor_forward_reverse_read_attributes[0]
    variant_position = tumor_forward_reverse_read_attributes[1]
    alt_allele_forward_reads = int(tumor_forward_reverse_read_attributes[2].split(',')[2])
    alt_allele_reverse_reads = int(tumor_forward_reverse_read_attributes[2].split(',')[3])

    if alt_allele_forward_reads < alt_allele_reverse_reads:
        forward_reverse_score = round(alt_allele_forward_reads / alt_allele_reverse_reads, ndigits=7)
    else:
        forward_reverse_score = round(alt_allele_reverse_reads / alt_allele_forward_reads, ndigits=7)

    print(chromosome, variant_position, forward_reverse_score, sep='\t')