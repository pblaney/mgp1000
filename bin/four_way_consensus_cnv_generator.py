#!/usr/bin/env python

"""
This script generates a consensus CNV BED file based on the per-segment agreement between the four
CNV callers included in the MGP1000 pipeline.
"""
import sys
import re

class CopyNumberSegment(object):

    def __init__(self, cn_segment_line):
        """Parse merged copy number segment file to construct object"""
        chrom, start, end, battenberg_cn, controlfreec_cn, sclust_cn, facets_cn = cn_segment_line.rstrip().split('\t')
        self.chrom = chrom
        self.start  = start
        self.end = end
        self.battenberg_cn = battenberg_cn
        self.controlfreec_cn = controlfreec_cn
        self.sclust_cn = sclust_cn
        self.facets_cn = facets_cn
        self.na_count = (battenberg_cn, controlfreec_cn, sclust_cn, facets_cn).count('.')
        self.segments_dict = {"battenberg": battenberg_cn,
                              "controlfreec": controlfreec_cn,
                              "sclust": sclust_cn,
                              "facets": facets_cn}

def consensus_writer(consensus_segment_tuple):
    """Write output BED line with consensus of per segment"""
    return '\t'.join(consensus_segment_tuple)

input_args = sys.argv

consensus_cnv_bed = open(input_args[2], 'w')

header = ("chrom", "start", "end",
          "consensus_cn", "caller_agreement",
          "battenberg_cn", "controlfreec_cn", "sclust_cn", "facets_cn")
consensus_cnv_bed.write('{0}\n'.format(consensus_writer(header)))

with open(input_args[1]) as merged_copy_number_file:
    for line in merged_copy_number_file:
        copy_num_obj = CopyNumberSegment(line)

        # First, check for how many tools generated a call for the segment
        if copy_num_obj.na_count == 3:

            # Output total copy number is equal to call from single tool
            single_caller_cn = [(key,value) for (key,value) in copy_num_obj.segments_dict.items() if value != '.'][0]
            single_caller_cn_tuple = (copy_num_obj.chrom, copy_num_obj.start, copy_num_obj.end,
                                      single_caller_cn[1], single_caller_cn[0],
                                      copy_num_obj.battenberg_cn, copy_num_obj.controlfreec_cn, copy_num_obj.sclust_cn, copy_num_obj.facets_cn)
            consensus_cnv_bed.write('{0}\n'.format(consensus_writer(single_caller_cn_tuple)))

        elif copy_num_obj.na_count == 2:
            double_caller_cn = [(key,value) for (key,value) in copy_num_obj.segments_dict.items() if value != '.']

            # If more than 1 tool generated a call for the segment, determine if the calls match
            if double_caller_cn[0][1] == double_caller_cn[1][1]:
                double_caller_agreement_cn_callers = ','.join([double_caller_cn[0][0], double_caller_cn[1][0]])
                double_caller_agreement_cn_tuple = (copy_num_obj.chrom, copy_num_obj.start, copy_num_obj.end,
                                                    double_caller_cn[0][1], double_caller_agreement_cn_callers,
                                                    copy_num_obj.battenberg_cn, copy_num_obj.controlfreec_cn, copy_num_obj.sclust_cn, copy_num_obj.facets_cn)
                consensus_cnv_bed.write('{0}\n'.format(consensus_writer(double_caller_agreement_cn_tuple)))

            else:
                double_caller_no_agreement_cn_tuple = (copy_num_obj.chrom, copy_num_obj.start, copy_num_obj.end,
                                                       ".", "no_agreement",
                                                       copy_num_obj.battenberg_cn, copy_num_obj.controlfreec_cn, copy_num_obj.sclust_cn, copy_num_obj.facets_cn)
                consensus_cnv_bed.write('{0}\n'.format(consensus_writer(double_caller_no_agreement_cn_tuple)))

        # If 3 tools generated a call for the segment, determine if the calls match
        elif copy_num_obj.na_count == 1:
            triple_caller_cn_dict = {}

            for caller,cn in copy_num_obj.segments_dict.items():
                if cn != ".":
                    if cn not in triple_caller_cn_dict:
                        triple_caller_cn_dict[cn] = [caller]
                    else:
                        triple_caller_cn_dict[cn].append(caller)

            # Catch segments with full agreement between 3 tools
            if len(triple_caller_cn_dict.keys()) == 1:
                triple_caller_full_agreement_cn_callers = ','.join(*triple_caller_cn_dict.values())
                triple_caller_full_agreement_cn_tuple = (copy_num_obj.chrom, copy_num_obj.start, copy_num_obj.end,
                                                         *triple_caller_cn_dict, triple_caller_full_agreement_cn_callers,
                                                         copy_num_obj.battenberg_cn, copy_num_obj.controlfreec_cn, copy_num_obj.sclust_cn, copy_num_obj.facets_cn)
                consensus_cnv_bed.write('{0}\n'.format(consensus_writer(triple_caller_full_agreement_cn_tuple)))

            # Catch segments with split agreement where 2/3 tools agree
            elif len(triple_caller_cn_dict.keys()) == 2:
                triple_partial_agreement_cn = [*triple_caller_cn_dict.keys()]
                triple_partial_agreement_cn_callers = [*triple_caller_cn_dict.values()]

                if len(triple_partial_agreement_cn_callers[0]) > len(triple_partial_agreement_cn_callers[1]):
                    triple_caller_partial_agreement_pair1_cn_callers = ','.join(triple_partial_agreement_cn_callers[0])
                    triple_caller_partial_agreement_pair1_cn_tuple = (copy_num_obj.chrom, copy_num_obj.start, copy_num_obj.end,
                                                                      triple_partial_agreement_cn[0], triple_caller_partial_agreement_pair1_cn_callers,
                                                                      copy_num_obj.battenberg_cn, copy_num_obj.controlfreec_cn, copy_num_obj.sclust_cn, copy_num_obj.facets_cn)
                    consensus_cnv_bed.write('{0}\n'.format(consensus_writer(triple_caller_partial_agreement_pair1_cn_tuple)))

                else:
                    triple_caller_partial_agreement_pair2_cn_callers = ','.join(triple_partial_agreement_cn_callers[1])
                    triple_caller_partial_agreement_pair2_cn_tuple = (copy_num_obj.chrom, copy_num_obj.start, copy_num_obj.end,
                                                                      triple_partial_agreement_cn[1], triple_caller_partial_agreement_pair2_cn_callers,
                                                                      copy_num_obj.battenberg_cn, copy_num_obj.controlfreec_cn, copy_num_obj.sclust_cn, copy_num_obj.facets_cn)
                    consensus_cnv_bed.write('{0}\n'.format(consensus_writer(triple_caller_partial_agreement_pair2_cn_tuple)))

            # Catch segment with no agreement between the 3 tools
            else:
                triple_caller_no_agreement_cn_tuple = (copy_num_obj.chrom, copy_num_obj.start, copy_num_obj.end,
                                                       ".", "no_agreement",
                                                       copy_num_obj.battenberg_cn, copy_num_obj.controlfreec_cn, copy_num_obj.sclust_cn, copy_num_obj.facets_cn)
                consensus_cnv_bed.write('{0}\n'.format(consensus_writer(triple_caller_no_agreement_cn_tuple)))

        # If 4 tools generated a call for the segment, determine if the calls match
        elif copy_num_obj.na_count == 0:
            quad_caller_cn_dict = {}

            for caller,cn in copy_num_obj.segments_dict.items():
                if cn not in quad_caller_cn_dict:
                    quad_caller_cn_dict[cn] = [caller]
                else:
                    quad_caller_cn_dict[cn].append(caller)

            # Catch segments with full agreement between 4 tools
            if len(quad_caller_cn_dict.keys()) == 1:
                quad_caller_full_agreement_cn_callers = ','.join(*quad_caller_cn_dict.values())
                quad_caller_full_agreement_cn_tuple = (copy_num_obj.chrom, copy_num_obj.start, copy_num_obj.end,
                                                       *quad_caller_cn_dict, quad_caller_full_agreement_cn_callers,
                                                       copy_num_obj.battenberg_cn, copy_num_obj.controlfreec_cn, copy_num_obj.sclust_cn, copy_num_obj.facets_cn)
                consensus_cnv_bed.write('{0}\n'.format(consensus_writer(quad_caller_full_agreement_cn_tuple)))

            # Catch segments with split agreement
            elif len(quad_caller_cn_dict.keys()) == 2:
                quad_partial_agreement_cn = [*quad_caller_cn_dict.keys()]
                quad_partial_agreement_cn_callers = [*quad_caller_cn_dict.values()]

                if len(quad_partial_agreement_cn_callers[0]) > len(quad_partial_agreement_cn_callers[1]):
                    quad_partial_agreement_pair1_cn_callers = ','.join(quad_partial_agreement_cn_callers[0])
                    quad_partial_agreement_pair1_cn_tuple = (copy_num_obj.chrom, copy_num_obj.start, copy_num_obj.end,
                                                             quad_partial_agreement_cn[0], quad_partial_agreement_pair1_cn_callers,
                                                             copy_num_obj.battenberg_cn, copy_num_obj.controlfreec_cn, copy_num_obj.sclust_cn, copy_num_obj.facets_cn)
                    consensus_cnv_bed.write('{0}\n'.format(consensus_writer(quad_partial_agreement_pair1_cn_tuple)))

                elif len(quad_partial_agreement_cn_callers[0]) < len(quad_partial_agreement_cn_callers[1]):
                    quad_partial_agreement_pair2_cn_callers = ','.join(quad_partial_agreement_cn_callers[1])
                    quad_partial_agreement_pair2_cn_tuple = (copy_num_obj.chrom, copy_num_obj.start, copy_num_obj.end,
                                                             quad_partial_agreement_cn[1], quad_partial_agreement_pair2_cn_callers,
                                                             copy_num_obj.battenberg_cn, copy_num_obj.controlfreec_cn, copy_num_obj.sclust_cn, copy_num_obj.facets_cn)
                    consensus_cnv_bed.write('{0}\n'.format(consensus_writer(quad_partial_agreement_pair2_cn_tuple)))

                # For segments with split 2 vs 2 agreement, prioritize the set with Battenberg
                else:
                    if re.search(r"battenberg", ','.join(quad_partial_agreement_cn_callers[0])):
                        quad_even_split_agreement_pair1_cn_callers = ','.join(quad_partial_agreement_cn_callers[0])
                        quad_even_split_agreement_pair1_cn_tuple = (copy_num_obj.chrom, copy_num_obj.start, copy_num_obj.end,
                                                                    quad_partial_agreement_cn[0], quad_even_split_agreement_pair1_cn_callers,
                                                                    copy_num_obj.battenberg_cn, copy_num_obj.controlfreec_cn, copy_num_obj.sclust_cn, copy_num_obj.facets_cn)
                        consensus_cnv_bed.write('{0}\n'.format(consensus_writer(quad_even_split_agreement_pair1_cn_tuple)))

                    elif re.search(r"battenberg", ','.join(quad_partial_agreement_cn_callers[1])):
                        quad_even_split_agreement_pair2_cn_callers = ','.join(quad_partial_agreement_cn_callers[1])
                        quad_even_split_agreement_pair2_cn_tuple = (copy_num_obj.chrom, copy_num_obj.start, copy_num_obj.end,
                                                                    quad_partial_agreement_cn[1], quad_even_split_agreement_pair2_cn_callers,
                                                                    copy_num_obj.battenberg_cn, copy_num_obj.controlfreec_cn, copy_num_obj.sclust_cn, copy_num_obj.facets_cn)
                        consensus_cnv_bed.write('{0}\n'.format(consensus_writer(quad_even_split_agreement_pair2_cn_tuple)))

            # Catch segments with lesser split agreement
            elif len(quad_caller_cn_dict.keys()) == 3:
                quad_lesser_partial_agreement_cn = [*quad_caller_cn_dict.keys()]
                quad_lesser_partial_agreement_cn_callers = [*quad_caller_cn_dict.values()]

                if len(quad_lesser_partial_agreement_cn_callers[0]) > len(quad_lesser_partial_agreement_cn_callers[1]) and len(quad_lesser_partial_agreement_cn_callers[0]) > len(quad_lesser_partial_agreement_cn_callers[2]):
                    quad_lesser_partial_agreement_pair1_cn_callers = ','.join(quad_lesser_partial_agreement_cn_callers[0])
                    quad_lesser_partial_agreement_pair1_cn_tuple = (copy_num_obj.chrom, copy_num_obj.start, copy_num_obj.end,
                                                                    quad_lesser_partial_agreement_cn[0], quad_lesser_partial_agreement_pair1_cn_callers,
                                                                    copy_num_obj.battenberg_cn, copy_num_obj.controlfreec_cn, copy_num_obj.sclust_cn, copy_num_obj.facets_cn)
                    consensus_cnv_bed.write('{0}\n'.format(consensus_writer(quad_lesser_partial_agreement_pair1_cn_tuple)))

                elif len(quad_lesser_partial_agreement_cn_callers[1]) > len(quad_lesser_partial_agreement_cn_callers[0]) and len(quad_lesser_partial_agreement_cn_callers[1]) > len(quad_lesser_partial_agreement_cn_callers[2]):
                    quad_lesser_partial_agreement_pair2_cn_callers = ','.join(quad_lesser_partial_agreement_cn_callers[1])
                    quad_lesser_partial_agreement_pair2_cn_tuple = (copy_num_obj.chrom, copy_num_obj.start, copy_num_obj.end,
                                                                    quad_lesser_partial_agreement_cn[1], quad_lesser_partial_agreement_pair2_cn_callers,
                                                                    copy_num_obj.battenberg_cn, copy_num_obj.controlfreec_cn, copy_num_obj.sclust_cn, copy_num_obj.facets_cn)
                    consensus_cnv_bed.write('{0}\n'.format(consensus_writer(quad_lesser_partial_agreement_pair2_cn_tuple)))

                elif len(quad_lesser_partial_agreement_cn_callers[2]) > len(quad_lesser_partial_agreement_cn_callers[0]) and len(quad_lesser_partial_agreement_cn_callers[2]) > len(quad_lesser_partial_agreement_cn_callers[1]):
                    quad_lesser_partial_agreement_pair3_cn_callers = ','.join(quad_lesser_partial_agreement_cn_callers[2])
                    quad_lesser_partial_agreement_pair3_cn_tuple = (copy_num_obj.chrom, copy_num_obj.start, copy_num_obj.end,
                                                                    quad_lesser_partial_agreement_cn[2], quad_lesser_partial_agreement_pair3_cn_callers,
                                                                    copy_num_obj.battenberg_cn, copy_num_obj.controlfreec_cn, copy_num_obj.sclust_cn, copy_num_obj.facets_cn)
                    consensus_cnv_bed.write('{0}\n'.format(consensus_writer(quad_lesser_partial_agreement_pair3_cn_tuple)))

            # Catch segments with no agreement between all 4 tools
            elif len(quad_caller_cn_dict.keys()) == 4:
                quad_caller_no_agreement_cn_tuple = (copy_num_obj.chrom, copy_num_obj.start, copy_num_obj.end,
                                                     ".", "no_agreement",
                                                     copy_num_obj.battenberg_cn, copy_num_obj.controlfreec_cn, copy_num_obj.sclust_cn, copy_num_obj.facets_cn)
                consensus_cnv_bed.write('{0}\n'.format(consensus_writer(quad_caller_no_agreement_cn_tuple)))

consensus_cnv_bed.close()
