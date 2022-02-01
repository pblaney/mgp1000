#!/usr/bin/env python

"""
This script generates a per segment BED file for total copy number and major/minor allele derived from
Control-FREEC's CNV/allele merged output.
There is a major assumption implemented here: for any segment with unassigned major/minor allele status,
a conservative estimation of heterozygosity is used.
"""
import sys
import math

class ControlFreecSegment(object):

    def __init__(self, control_freec_line):
        """Parse Control-FREEC allele per segment file to construct object"""
        chrom, start, end, total_cn, genotype = control_freec_line.rstrip().split('\t')
        self.chrom = chrom
        self.start = start
        self.end = end
        self.total_cn = total_cn
        self.genotype = genotype

def output_segment_writer(segment_tuple):
    """Write output BED line with total copy number and alleles per segment"""
    return '\t'.join(segment_tuple)

input_args = sys.argv

control_freec_alleles_bed = open(input_args[2], 'w')
control_freec_cnv_bed = open(input_args[3], 'w')

with open(input_args[1]) as control_freec_mapped_alleles:
    for line in control_freec_mapped_alleles:
        control_freec_segment_obj = ControlFreecSegment(line)

        cnv_tup = (control_freec_segment_obj.chrom, control_freec_segment_obj.start, control_freec_segment_obj.end, control_freec_segment_obj.total_cn)
        control_freec_cnv_bed.write('{0}\n'.format(output_segment_writer(cnv_tup)))

        # First fix all unassigned genotype segments
        if control_freec_segment_obj.genotype == "-" or control_freec_segment_obj.genotype == ".":

            if control_freec_segment_obj.total_cn == "0":
                bialleleic_deletion_inferred_allele_tup = (control_freec_segment_obj.chrom, control_freec_segment_obj.start, control_freec_segment_obj.end, "0/0")
                control_freec_alleles_bed.write('{0}\n'.format(output_segment_writer(bialleleic_deletion_inferred_allele_tup)))

            elif control_freec_segment_obj.total_cn == "1":
                monoalleleic_deletion_inferred_allele_tup = (control_freec_segment_obj.chrom, control_freec_segment_obj.start, control_freec_segment_obj.end, "1/0")
                control_freec_alleles_bed.write('{0}\n'.format(output_segment_writer(monoalleleic_deletion_inferred_allele_tup)))

            else:
                # Use known total copy number value to determine the heterozygous major/minor allele assignment
                cn_per_allele = int(int(control_freec_segment_obj.total_cn) / 2)
                is_unbalanced = round(math.fmod(int(control_freec_segment_obj.total_cn), 2))

                if is_unbalanced == 0:
                    balanced_heterozygous_inferred_alleles = '/'.join((str(cn_per_allele), str(cn_per_allele)))
                    balanced_heterozygous_inferred_allele_tup = (control_freec_segment_obj.chrom, control_freec_segment_obj.start, control_freec_segment_obj.end, balanced_heterozygous_inferred_alleles)
                    control_freec_alleles_bed.write('{0}\n'.format(output_segment_writer(balanced_heterozygous_inferred_allele_tup)))

                elif is_unbalanced == 1:
                    unbalanced_heterozygous_inferred_alleles = '/'.join((str(cn_per_allele + is_unbalanced), str(cn_per_allele)))
                    unbalanced_heterozygous_inferred_allele_tup = (control_freec_segment_obj.chrom, control_freec_segment_obj.start, control_freec_segment_obj.end, unbalanced_heterozygous_inferred_alleles)
                    control_freec_alleles_bed.write('{0}\n'.format(output_segment_writer(unbalanced_heterozygous_inferred_allele_tup)))

        else:
            # Convert the 'AB' convention to numeric '1/1' allele format
            major_minor_alleles = '/'.join([str(control_freec_segment_obj.genotype.count("A")), str(control_freec_segment_obj.genotype.count("B"))])
            allele_segment_tup = (control_freec_segment_obj.chrom, control_freec_segment_obj.start, control_freec_segment_obj.end, major_minor_alleles)
            control_freec_alleles_bed.write('{0}\n'.format(output_segment_writer(allele_segment_tup)))

control_freec_alleles_bed.close()
control_freec_cnv_bed.close()
