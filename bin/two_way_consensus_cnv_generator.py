#!/usr/bin/env python

"""
This script generates a consensus CNV BED file based on the per-segment agreement between the two
CNV callers included in the MGP1000 pipeline.
"""
import sys

class CopyNumberSegment(object):

    def __init__(self, cn_segment_line):
        """Parse merged copy number segment file to construct object"""
        chrom, start, end, battenberg_cn, facets_cn = cn_segment_line.rstrip().split('\t')
        self.chrom = chrom
        self.start  = start
        self.end = end
        self.battenberg_cn = battenberg_cn
        self.facets_cn = facets_cn
        self.na_count = (battenberg_cn, facets_cn).count('.')
        self.segments_dict = {"battenberg": battenberg_cn,
                              "facets": facets_cn}

def consensus_writer(consensus_segment_tuple):
    """Write output BED line with consensus of per segment"""
    return '\t'.join(consensus_segment_tuple)

input_args = sys.argv

consensus_cnv_bed = open(input_args[2], 'w')

header = ("chrom", "start", "end",
          "consensus_cn", "caller_agreement",
          "battenberg_cn", "facets_cn")
consensus_cnv_bed.write('{0}\n'.format(consensus_writer(header)))

with open(input_args[1]) as merged_copy_number_file:
    for line in merged_copy_number_file:
        copy_num_obj = CopyNumberSegment(line)

        # First, check for how many tools generated a call for the segment
        if copy_num_obj.na_count == 1:

            # Output total copy number is equal to call from single tool
            single_caller_cn = [(key,value) for (key,value) in copy_num_obj.segments_dict.items() if value != '.'][0]
            single_caller_cn_tuple = (copy_num_obj.chrom, copy_num_obj.start, copy_num_obj.end,
                                      single_caller_cn[1], single_caller_cn[0],
                                      copy_num_obj.battenberg_cn, copy_num_obj.facets_cn)
            consensus_cnv_bed.write('{0}\n'.format(consensus_writer(single_caller_cn_tuple)))


        elif copy_num_obj.na_count == 0:
            double_caller_cn = [(key,value) for (key,value) in copy_num_obj.segments_dict.items() if value != '.']

            # If more than 1 tool generated a call for the segment, determine if the calls match
            if double_caller_cn[0][1] == double_caller_cn[1][1]:
                double_caller_agreement_cn_callers = ','.join([double_caller_cn[0][0], double_caller_cn[1][0]])
                double_caller_agreement_cn_tuple = (copy_num_obj.chrom, copy_num_obj.start, copy_num_obj.end,
                                                    double_caller_cn[0][1], double_caller_agreement_cn_callers,
                                                    copy_num_obj.battenberg_cn, copy_num_obj.facets_cn)
                consensus_cnv_bed.write('{0}\n'.format(consensus_writer(double_caller_agreement_cn_tuple)))

            else:
                double_caller_no_agreement_cn_tuple = (copy_num_obj.chrom, copy_num_obj.start, copy_num_obj.end,
                                                       ".", "no_agreement",
                                                       copy_num_obj.battenberg_cn, copy_num_obj.facets_cn)
                consensus_cnv_bed.write('{0}\n'.format(consensus_writer(double_caller_no_agreement_cn_tuple)))

consensus_cnv_bed.close()
