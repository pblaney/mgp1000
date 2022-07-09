#!/usr/bin/env python

"""
This script extracts the clonal and subclonal copy number profiles from Battenberg's output
"""
import sys
import math

class BattenbergSegment(object):

    def __init__(self, battenberg_line):
        """Parse Battenberg file to construct object"""
        chrom, start, end, state1_maj_allele, state1_min_allele, state1_fraction, state2_maj_allele, state2_min_allele, state2_fraction = battenberg_line.rstrip().split('\t')
        self.chrom = chrom
        self.start = start
        self.end = end
        self.state1_maj_allele = state1_maj_allele
        self.state1_min_allele = state1_min_allele
        self.state1_fraction = float(state1_fraction)
        self.state2_maj_allele = state2_maj_allele
        self.state2_min_allele = state2_min_allele
        self.state2_fraction = float(state2_fraction) if state2_fraction != "NA" else "NA"

def output_writer(segment_tuple):
    """Write output to designated file"""
    return '\t'.join(segment_tuple)

input_args = sys.argv

clonal_cn_output = open(input_args[2], 'w')

clonal_allele_output = open(input_args[3], 'w')

subclonal_output = open(input_args[4], 'w')
subclonal_header = ("chrom", "start", "end",
                    "subclonal_major_allele", "subclonal_minor_allele",
                    "subclonal_fraction")
subclonal_output.write('{0}\n'.format(output_writer(subclonal_header)))

with open(input_args[1]) as battenberg_cn_profile:
    for line in battenberg_cn_profile:
        battenberg_segment_obj = BattenbergSegment(line)

        # Identify the clonal output determined by the highest state fraction of the two, then separately write output to
        # total CN and allele files for consensus
        if battenberg_segment_obj.state2_fraction == "NA":
            full_clonal_total_cn = str(sum([int(battenberg_segment_obj.state1_maj_allele), int(battenberg_segment_obj.state1_min_allele)]))
            full_clonal_cn_tuple = (battenberg_segment_obj.chrom, battenberg_segment_obj.start, battenberg_segment_obj.end,
                                    full_clonal_total_cn)
            clonal_cn_output.write('{0}\n'.format(output_writer(full_clonal_cn_tuple)))

            full_clonal_allele_tuple = (battenberg_segment_obj.chrom, battenberg_segment_obj.start, battenberg_segment_obj.end,
                                        '/'.join([battenberg_segment_obj.state1_maj_allele, battenberg_segment_obj.state1_min_allele]))
            clonal_allele_output.write('{0}\n'.format(output_writer(full_clonal_allele_tuple)))

        elif battenberg_segment_obj.state1_fraction > battenberg_segment_obj.state2_fraction:
            state1_clonal_total_cn = str(sum([int(battenberg_segment_obj.state1_maj_allele), int(battenberg_segment_obj.state1_min_allele)]))
            state1_clonal_cn_tuple = (battenberg_segment_obj.chrom, battenberg_segment_obj.start, battenberg_segment_obj.end,
                                      state1_clonal_total_cn)
            clonal_cn_output.write('{0}\n'.format(output_writer(state1_clonal_cn_tuple)))

            state1_clonal_allele_tuple = (battenberg_segment_obj.chrom, battenberg_segment_obj.start, battenberg_segment_obj.end,
                                          '/'.join([battenberg_segment_obj.state1_maj_allele, battenberg_segment_obj.state1_min_allele]))
            clonal_allele_output.write('{0}\n'.format(output_writer(state1_clonal_allele_tuple)))

            # Write subclonal segment to separate input
            state2_subclonal_tuple = (battenberg_segment_obj.chrom, battenberg_segment_obj.start, battenberg_segment_obj.end,
                                      battenberg_segment_obj.state2_maj_allele, battenberg_segment_obj.state2_min_allele, str(battenberg_segment_obj.state2_fraction))
            subclonal_output.write('{0}\n'.format(output_writer(state2_subclonal_tuple)))

        elif battenberg_segment_obj.state2_fraction > battenberg_segment_obj.state1_fraction:
            state2_clonal_total_cn = str(sum([int(battenberg_segment_obj.state2_maj_allele), int(battenberg_segment_obj.state2_min_allele)]))
            state2_clonal_cn_tuple = (battenberg_segment_obj.chrom, battenberg_segment_obj.start, battenberg_segment_obj.end,
                                      state2_clonal_total_cn)
            clonal_cn_output.write('{0}\n'.format(output_writer(state2_clonal_cn_tuple)))

            state2_clonal_allele_tuple = (battenberg_segment_obj.chrom, battenberg_segment_obj.start, battenberg_segment_obj.end,
                                          '/'.join([battenberg_segment_obj.state2_maj_allele, battenberg_segment_obj.state2_min_allele]))
            clonal_allele_output.write('{0}\n'.format(output_writer(state2_clonal_allele_tuple)))

            # Write subclonal segment to separate input
            state1_subclonal_tuple = (battenberg_segment_obj.chrom, battenberg_segment_obj.start, battenberg_segment_obj.end,
                                      battenberg_segment_obj.state1_maj_allele, battenberg_segment_obj.state2_min_allele, str(battenberg_segment_obj.state1_fraction))
            subclonal_output.write('{0}\n'.format(output_writer(state1_subclonal_tuple)))
        else:
            continue

clonal_cn_output.close()
clonal_allele_output.close()
subclonal_output.close()
