#!/usr/bin/env python

"""
This script generates a consensus CNV BED file based on the per-segment agreement between the three
CNV callers included in the MGP1000 pipeline.
"""
import sys

class CopyNumberSegment(object):

    def __init__(self, cn_segment_line):
        """Parse merged copy number segment file to construct object"""
        chrom, start, end, ascat_cn, controlfreec_cn, sclust_cn = cn_segment_line.rstrip().split('\t')
        self.chrom = chrom
        self.start  = start
        self.end = end
        self.ascat_cn = ascat_cn
        self.controlfreec_cn = controlfreec_cn
        self.sclust_cn = sclust_cn
        self.na_count = (ascat_cn, controlfreec_cn, sclust_cn).count('NA')
        self.segments_dict = {"ascat": ascat_cn,
                              "controlfreec": controlfreec_cn,
                              "sclust": sclust_cn}

def consensus_writer(consensus_cnv_segment_tuple):
    '''Return output BED line with CNV consensus of per segment'''
    return '\t'.join(consensus_cnv_segment_tuple)

input_args = sys.argv

consensus_cnv_bed = open(input_args[2], 'w')

header = ("chrom", "start", "end",
          "consensus_total_cn", "caller_agreement_cn",
          "ascat_total_cn", "controlfreec_total_cn", "sclust_total_cn")
consensus_cnv_bed.write('{0}\n'.format(consensus_writer(header)))

with open(input_args[1]) as merged_copy_number_file:
    for line in merged_copy_number_file:
        copy_num_obj = CopyNumberSegment(line)

        if copy_num_obj.na_count == 2:
            single_caller_cn = [(key,value) for (key,value) in copy_num_obj.segments_dict.items() if value != 'NA'][0]
            single_caller_cn_tuple = (copy_num_obj.chrom, copy_num_obj.start, copy_num_obj.end,
                                      single_caller_cn[1], single_caller_cn[0],
                                      copy_num_obj.ascat_cn, copy_num_obj.controlfreec_cn, copy_num_obj.sclust_cn)
            consensus_cnv_bed.write('{0}\n'.format(consensus_writer(single_caller_cn_tuple)))

        elif copy_num_obj.na_count == 1:
            double_caller_cn = [(key,value) for (key,value) in copy_num_obj.segments_dict.items() if value != 'NA']

            if double_caller_cn[0][1] == double_caller_cn[1][1]:
                double_caller_cn_callers = ','.join([double_caller_cn[0][0], double_caller_cn[1][0]])
                double_caller_cn_tuple = (copy_num_obj.chrom, copy_num_obj.start, copy_num_obj.end,
                                          double_caller_cn[0][1], double_caller_cn_callers,
                                          copy_num_obj.ascat_cn, copy_num_obj.controlfreec_cn, copy_num_obj.sclust_cn)
                consensus_cnv_bed.write('{0}\n'.format(consensus_writer(double_caller_cn_tuple)))
            else:
                double_caller_no_agreement_cn_tuple = (copy_num_obj.chrom, copy_num_obj.start, copy_num_obj.end,
                                                       "-", "no_agreement",
                                                       copy_num_obj.ascat_cn, copy_num_obj.controlfreec_cn, copy_num_obj.sclust_cn)
                consensus_cnv_bed.write('{0}\n'.format(consensus_writer(double_caller_no_agreement_cn_tuple)))

        elif copy_num_obj.na_count == 0:
            triple_caller_cn_dict = {}

            for caller,cn in copy_num_obj.segments_dict.items():
                if cn not in triple_caller_cn_dict:
                    triple_caller_cn_dict[cn] = [caller]
                else:
                    triple_caller_cn_dict[cn].append(caller)

            if len(triple_caller_cn_dict.keys()) == 1:
                triple_caller_cn_callers = ','.join(*triple_caller_cn_dict.values())
                triple_caller_cn_tuple = (copy_num_obj.chrom, copy_num_obj.start, copy_num_obj.end,
                                          *triple_caller_cn_dict, triple_caller_cn_callers,
                                          copy_num_obj.ascat_cn, copy_num_obj.controlfreec_cn, copy_num_obj.sclust_cn)
                consensus_cnv_bed.write('{0}\n'.format(consensus_writer(triple_caller_cn_tuple)))

            elif len(triple_caller_cn_dict.keys()) == 2:
                triple_split_cn = [*triple_caller_cn_dict.keys()]
                triple_split_cn_callers = [*triple_caller_cn_dict.values()]

                if len(triple_split_cn_callers[0]) > len(triple_split_cn_callers[1]):
                    pair1_of_three_cn_callers = ','.join(triple_split_cn_callers[0])
                    pair1_of_three_cn_tuple = (copy_num_obj.chrom, copy_num_obj.start, copy_num_obj.end,
                                               triple_split_cn[0], pair1_of_three_cn_callers,
                                               copy_num_obj.ascat_cn, copy_num_obj.controlfreec_cn, copy_num_obj.sclust_cn)
                    consensus_cnv_bed.write('{0}\n'.format(consensus_writer(pair1_of_three_cn_tuple)))
                else:
                    pair2_of_three_cn_callers = ','.join(triple_split_cn_callers[1])
                    pair2_of_three_cn_tuple = (copy_num_obj.chrom, copy_num_obj.start, copy_num_obj.end,
                                               triple_split_cn[1], pair2_of_three_cn_callers,
                                               copy_num_obj.ascat_cn, copy_num_obj.controlfreec_cn, copy_num_obj.sclust_cn)
                    consensus_cnv_bed.write('{0}\n'.format(consensus_writer(pair2_of_three_cn_tuple)))

            else:
                triple_caller_no_agreement_cn_tuple = (copy_num_obj.chrom, copy_num_obj.start, copy_num_obj.end,
                                                       "-", "no_agreement",
                                                       copy_num_obj.ascat_cn, copy_num_obj.controlfreec_cn, copy_num_obj.sclust_cn)
                consensus_cnv_bed.write('{0}\n'.format(consensus_writer(triple_caller_no_agreement_cn_tuple)))

consensus_cnv_bed.close()
