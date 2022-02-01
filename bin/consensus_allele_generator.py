#!/usr/bin/env python

"""
This script generates a consensus major/minor allele BED file based on the per-segment agreement between the three
CNV callers included in the MGP1000 pipeline.
"""
import sys

class AlleleSegment(object):

    def __init__(self, allele_segment_line):
        """Parse merged allele segment file to construct object"""
        chrom, start, end, ascat_alleles, controlfreec_alleles, sclust_alleles = allele_segment_line.rstrip().split('\t')
        self.chrom = chrom
        self.start  = start
        self.end = end
        self.ascat_alleles = ascat_alleles
        self.controlfreec_alleles = controlfreec_alleles
        self.sclust_alleles = sclust_alleles
        self.na_count = (ascat_alleles, controlfreec_alleles, sclust_alleles).count('.')
        self.alleles_dict = {"ascat": ascat_alleles,
                             "controlfreec": controlfreec_alleles,
                             "sclust": sclust_alleles}

def consensus_writer(consensus_allele_segment_tuple):
    """Return output BED line with allele consensus of per segment"""
    return '\t'.join(consensus_allele_segment_tuple)

input_args = sys.argv

consensus_allele_bed = open(input_args[2], 'w')

header = ("chrom", "start", "end",
          "consensus_major_allele", "consensus_minor_allele", "allele_caller_agreement",
          "ascat_alleles", "controlfreec_alleles", "sclust_alleles")
consensus_allele_bed.write('{0}\n'.format(consensus_writer(header)))

with open(input_args[1]) as merged_allele_file:
    for line in merged_allele_file:
        allele_obj = AlleleSegment(line)

        if allele_obj.na_count == 2:
            single_caller_allele = [(key,value) for (key,value) in allele_obj.alleles_dict.items() if value != '.'][0]
            single_caller_allele_tuple = (allele_obj.chrom, allele_obj.start, allele_obj.end,
                                          single_caller_allele[1].replace("/", "\t"), single_caller_allele[0],
                                          allele_obj.ascat_alleles, allele_obj.controlfreec_alleles, allele_obj.sclust_alleles)
            consensus_allele_bed.write('{0}\n'.format(consensus_writer(single_caller_allele_tuple)))

        elif allele_obj.na_count == 1:
            double_caller_allele = [(key,value) for (key,value) in allele_obj.alleles_dict.items() if value != '.']

            if double_caller_allele[0][1] == double_caller_allele[1][1]:
                double_caller_allele_callers = ','.join([double_caller_allele[0][0], double_caller_allele[1][0]])
                double_caller_allele_tuple = (allele_obj.chrom, allele_obj.start, allele_obj.end,
                                              double_caller_allele[0][1].replace("/", "\t"), double_caller_allele_callers,
                                              allele_obj.ascat_alleles, allele_obj.controlfreec_alleles, allele_obj.sclust_alleles)
                consensus_allele_bed.write('{0}\n'.format(consensus_writer(double_caller_allele_tuple)))
            else:
                double_caller_no_agreement_allele_tuple = (allele_obj.chrom, allele_obj.start, allele_obj.end,
                                                           "-", "-", "no_agreement",
                                                           allele_obj.ascat_alleles, allele_obj.controlfreec_alleles, allele_obj.sclust_alleles)
                consensus_allele_bed.write('{0}\n'.format(consensus_writer(double_caller_no_agreement_allele_tuple)))

        elif allele_obj.na_count == 0:
            triple_caller_allele_dict = {}

            for caller,allele in allele_obj.alleles_dict.items():
                if allele not in triple_caller_allele_dict:
                    triple_caller_allele_dict[allele] = [caller]
                else:
                    triple_caller_allele_dict[allele].append(caller)

            if len(triple_caller_allele_dict.keys()) == 1:
                triple_caller_allele_callers = ','.join(*triple_caller_allele_dict.values())
                triple_caller_allele_tuple = (allele_obj.chrom, allele_obj.start, allele_obj.end,
                                              '{0}'.format(*triple_caller_allele_dict).replace("/", "\t"), triple_caller_allele_callers,
                                              allele_obj.ascat_alleles, allele_obj.controlfreec_alleles, allele_obj.sclust_alleles)
                consensus_allele_bed.write('{0}\n'.format(consensus_writer(triple_caller_allele_tuple)))

            elif len(triple_caller_allele_dict.keys()) == 2:
                triple_split_allele = [*triple_caller_allele_dict.keys()]
                triple_split_allele_callers = [*triple_caller_allele_dict.values()]

                if len(triple_split_allele_callers[0]) > len(triple_split_allele_callers[1]):
                    pair1_of_three_allele_callers = ','.join(triple_split_allele_callers[0])
                    pair1_of_three_allele_tuple = (allele_obj.chrom, allele_obj.start, allele_obj.end,
                                                   triple_split_allele[0].replace("/", "\t"), pair1_of_three_allele_callers,
                                                   allele_obj.ascat_alleles, allele_obj.controlfreec_alleles, allele_obj.sclust_alleles)
                    consensus_allele_bed.write('{0}\n'.format(consensus_writer(pair1_of_three_allele_tuple)))
                else:
                    pair2_of_three_allele_callers = ','.join(triple_split_allele_callers[1])
                    pair2_of_three_allele_tuple = (allele_obj.chrom, allele_obj.start, allele_obj.end,
                                                   triple_split_allele[1].replace("/", "\t"), pair2_of_three_allele_callers,
                                                   allele_obj.ascat_alleles, allele_obj.controlfreec_alleles, allele_obj.sclust_alleles)
                    consensus_allele_bed.write('{0}\n'.format(consensus_writer(pair2_of_three_allele_tuple)))
            else:
                triple_caller_no_agreement_allele_tuple = (allele_obj.chrom, allele_obj.start, allele_obj.end,
                                                           "-", "-", "no_agreement",
                                                           allele_obj.ascat_alleles, allele_obj.controlfreec_alleles, allele_obj.sclust_alleles)
                consensus_allele_bed.write('{0}\n'.format(consensus_writer(triple_caller_no_agreement_allele_tuple)))

consensus_allele_bed.close()
