#!/usr/bin/env python

"""
This script generates a BED file based the per CNV segment allele states derived from Control-FREEC's output files
"""
import sys

class ControlFreecAlleleSegment(object):

    def __init__(self, control_freec_line):
        """Parse Control-FREEC allele per segment file to construct object"""
        chrom, start, end, total_cn, genotype = control_freec_line.rstrip().split('\t')
        self.chrom = chrom
        self.start = start
        self.end = end
        self.total_cn = total_cn
        self.genotype = genotype

def allele_segment_writer(allele_segment_tuple):
    '''Write output BED line with alleles per segment'''
    return '\t'.join(allele_segment_tuple)

input_args = sys.argv

control_freec_alleles_bed = open(input_args[2], 'w')

with open(input_args[1]) as control_freec_mapped_alleles:
    for line in control_freec_mapped_alleles:
        allele_segment_obj = ControlFreecAlleleSegment(line)

        if allele_segment_obj.genotype == "-":

            if allele_segment_obj.total_cn == "0":
                bialleleic_deletion_tuple = (allele_segment_obj.chrom, allele_segment_obj.start, allele_segment_obj.end, "0/0")
                control_freec_alleles_bed.write('{0}\n'.format(allele_segment_writer(bialleleic_deletion_tuple)))
            else:
                nonzero_cn_with_undetermined_genotype_tuple = (allele_segment_obj.chrom, allele_segment_obj.start, allele_segment_obj.end, "NA")
                control_freec_alleles_bed.write('{0}\n'.format(allele_segment_writer(nonzero_cn_with_undetermined_genotype_tuple)))

        elif allele_segment_obj.genotype == ".":

            if allele_segment_obj.total_cn == "1":
                inferred_monoallelic_deletion_tuple = (allele_segment_obj.chrom, allele_segment_obj.start, allele_segment_obj.end, "1/0")
                control_freec_alleles_bed.write('{0}\n'.format(allele_segment_writer(inferred_monoallelic_deletion_tuple)))
            elif allele_segment_obj.total_cn == "2":
                neutral_heterozygous_tuple = (allele_segment_obj.chrom, allele_segment_obj.start, allele_segment_obj.end, "1/1")
                control_freec_alleles_bed.write('{0}\n'.format(allele_segment_writer(neutral_heterozygous_tuple)))

        else:

            if allele_segment_obj.total_cn == "1" and allele_segment_obj.genotype == "A":
                monoallelic_deletion_tuple = (allele_segment_obj.chrom, allele_segment_obj.start, allele_segment_obj.end, "1/0")
                control_freec_alleles_bed.write('{0}\n'.format(allele_segment_writer(monoallelic_deletion_tuple)))
            elif allele_segment_obj.total_cn == "2" and allele_segment_obj.genotype == "AA":
                neutral_homozygous_tuple = (allele_segment_obj.chrom, allele_segment_obj.start, allele_segment_obj.end, "2/0")
                control_freec_alleles_bed.write('{0}\n'.format(allele_segment_writer(neutral_homozygous_tuple)))
            else:
                nonneutral_heterozygous_alleles = '/'.join([str(allele_segment_obj.genotype.count("A")), str(allele_segment_obj.genotype.count("B"))])
                nonneutral_heterozygous_tuple = (allele_segment_obj.chrom, allele_segment_obj.start, allele_segment_obj.end, nonneutral_heterozygous_alleles)
                control_freec_alleles_bed.write('{0}\n'.format(allele_segment_writer(nonneutral_heterozygous_tuple)))
