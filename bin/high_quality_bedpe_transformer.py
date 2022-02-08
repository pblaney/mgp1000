#!/usr/bin/env python

"""
This script transforms an AnnotSV annotated set of consensus structural variants into a high quality BEDPE
"""
import sys

class AnnotStructuralVariant(object):

    def __init__(self, annot_sv_line):
        """Parse AnnotSV annotated BED file"""
        annot_sv_entry = annot_sv_line.rstrip().split('\t')
        self.chrom = annot_sv_entry[0]
        self.start = annot_sv_entry[1]
        self.end = annot_sv_entry[2]
        self.cytoband = annot_sv_entry[3]
        self.length = "Length=" + annot_sv_entry[4] if annot_sv_entry[4] != "0" else "Length=."
        self.type = annot_sv_entry[5]
        self.annotsv_id = "AnnotSV_id=" + annot_sv_entry[6]
        self.annotsv_score = "AnnotSV_score=" + annot_sv_entry[7]
        self.acmg_class = "ACMG_score=" + annot_sv_entry[8]
        self.vcf_id = annot_sv_entry[9]
        self.vcf_info = annot_sv_entry[10]
        self.delly_id = annot_sv_entry[11].split(':')[7] if annot_sv_entry[11].split(':')[7] != 'NaN' else '.'
        self.manta_id = annot_sv_entry[12].split(':')[7].replace('_',':') if annot_sv_entry[12].split(':')[7] != 'NaN' else '.'
        self.svaba_id = annot_sv_entry[13].split(':')[7].replace('_',':') if annot_sv_entry[13].split(':')[7] != 'NaN' else '.'
        self.genes = "Genes=" + annot_sv_entry[14].replace(';',',')
        self.gene_count = "Gene_count=" + annot_sv_entry[15]
        self.gc_content_end1 = "GC_content_end1=" + annot_sv_entry[16]
        self.gc_content_end2 = "GC_content_end2=" + annot_sv_entry[17]
        self.repeat_coord_end1 = "Repeat_coord_end1=" + annot_sv_entry[18].replace(';',',')
        self.repeat_type_end1 = "Repeat_type_end1=" + annot_sv_entry[19].replace(';',',')
        self.repeat_coord_end2 = "Repeat_coord_end2=" + annot_sv_entry[20].replace(';',',')
        self.repeat_type_end2 = "Repeat_type_end2=" + annot_sv_entry[21].replace(';',',')
        self.gap_end1 = "Gap_end1=" + annot_sv_entry[22].replace(';',',')
        self.gap_end2 = "Gap_end2=" + annot_sv_entry[23].replace(';',',')
        self.segdup_end1 = "SegDup_end1=" + annot_sv_entry[24].replace(';',',')
        self.segdup_end2 = "SegDup_end2=" + annot_sv_entry[25].replace(';',',')
        self.blacklist_end1 = "ENCODE_blacklist_end1=" + annot_sv_entry[26]
        self.blacklist_type_end1 = "ENCODE_blacklist_type_end1=" + annot_sv_entry[27].replace(' ','')
        self.blacklist_end2 = "ENCODE_blacklist_end2=" + annot_sv_entry[28]
        self.blacklist_type_end2 = "ENCODE_blacklist_type_end2=" + annot_sv_entry[29].replace(' ','')

        # Isolate caller vector field from INFO column and rewrite to be more informative
        callers_string = self.vcf_info.split(';')[1].replace('SUPP_VEC=','')
        caller_combo_dict = {'111': 'delly,manta,svaba', '100': 'delly', '010': 'manta', '001': 'svaba',
                             '110': 'delly,manta', '101': 'delly,svaba', '011': 'manta,svaba'}
        self.callers = "Callers=" + caller_combo_dict[callers_string]

        # Isolate details for breakend junction mate from INFO column
        self.partner_chrom = self.vcf_info.split(';')[5].replace('CHR2=','')
        self.partner_start = self.vcf_info.split(';')[6].replace('END=','')
        self.partner_end = str(int(self.partner_start) + 1)

        # Isolate strand information from INFO column
        strands_string = self.vcf_info.split(';')[9].replace('STRANDS=','')
        self.strand_1 = strands_string[0]
        self.strand_2 = strands_string[1]

def read_support_dictonary_creater(read_support_file):
    """Read in table of read support values from different callers"""
    read_support_dict = {}
    for read_supp_line in read_support_file:
        call_id, split_reads, discordant_reads = read_supp_line.rstrip().split('\t')
        read_support_dict[call_id] = [split_reads,discordant_reads]
    return read_support_dict

def bedpe_output_writer(bedpe_output_tuple):
    """Return BEDPE formatted output line"""
    return '\t'.join(bedpe_output_tuple)

input_args = sys.argv

# Create output file based on user-defined filename
high_quality_sv_bedpe = open(input_args[5], 'w')

column_header = ("#chrom1", "start1", "end1", "chrom2", "start2", "end2",
                 "type", "score", "strand1", "strand2",
                 "cytoband", "support", "annotations")
high_quality_sv_bedpe.write('{0}\n'.format(bedpe_output_writer(column_header)))

# Create dictionaries of read support data from each caller
with open(input_args[2]) as manta_read_support, open(input_args[3]) as svaba_read_support, open(input_args[4]) as delly_read_support:
    manta_read_support_dict = read_support_dictonary_creater(manta_read_support)
    svaba_read_support_dict = read_support_dictonary_creater(svaba_read_support)
    delly_read_support_dict = read_support_dictonary_creater(delly_read_support)

# Read in AnnotSV BED file and create BEDPE formatted entries
with open(input_args[1]) as annot_sv_bed:

    inter_chrom_list = []
    for sv_record in annot_sv_bed:
        annot_sv_obj = AnnotStructuralVariant(sv_record)

        # Create score column using combination of annotations
        score = ';'.join((annot_sv_obj.callers, annot_sv_obj.annotsv_score, annot_sv_obj.acmg_class, annot_sv_obj.length))

        # Create support column using read support information from each caller if available
        support_list = []
        if annot_sv_obj.delly_id != '.':
            delly_discordant_reads = "Delly_DR=" + delly_read_support_dict[annot_sv_obj.delly_id][0]
            delly_split_reads = "Delly_SR=" + delly_read_support_dict[annot_sv_obj.delly_id][1]
            support_list.append(','.join((delly_discordant_reads, delly_split_reads)))

        if annot_sv_obj.manta_id != '.':
            manta_discordant_reads = "Manta_DR=" + manta_read_support_dict[annot_sv_obj.manta_id][0]
            manta_split_reads = "Manta_SR=" + manta_read_support_dict[annot_sv_obj.manta_id][1]
            support_list.append(','.join((manta_discordant_reads, manta_split_reads)))

        if annot_sv_obj.svaba_id != '.':
            svaba_discordant_reads = "Svaba_DR=" + svaba_read_support_dict[annot_sv_obj.svaba_id][0]
            svaba_split_reads = "Svaba_SR=" + svaba_read_support_dict[annot_sv_obj.svaba_id][1]
            support_list.append(','.join((svaba_discordant_reads, svaba_split_reads)))

        support = ';'.join(tuple(support_list))

        # Create annotations column using combination of remaining annotations
        annotations = ';'.join((annot_sv_obj.genes, annot_sv_obj.gc_content_end1, annot_sv_obj.gc_content_end2,
                                annot_sv_obj.repeat_coord_end1, annot_sv_obj.repeat_type_end1,
                                annot_sv_obj.repeat_coord_end2, annot_sv_obj.repeat_type_end2,
                                annot_sv_obj.gap_end1, annot_sv_obj.gap_end2,
                                annot_sv_obj.segdup_end1, annot_sv_obj.segdup_end2,
                                annot_sv_obj.blacklist_end1, annot_sv_obj.blacklist_type_end1,
                                annot_sv_obj.blacklist_end2, annot_sv_obj.blacklist_type_end2,
                                annot_sv_obj.annotsv_id))

        # Main functionality, generate the BEDPE format output depending on intra- or inter- chromosomal SV type
        if annot_sv_obj.type == 'BND':

            bedpe_interchrom_entry_list = [annot_sv_obj.chrom, annot_sv_obj.start, annot_sv_obj.end,
                                           annot_sv_obj.partner_chrom, annot_sv_obj.partner_start, annot_sv_obj.partner_end,
                                           annot_sv_obj.type, score, annot_sv_obj.strand_1, annot_sv_obj.strand_2,
                                           annot_sv_obj.cytoband, support, annotations]
            inter_chrom_list.append(bedpe_interchrom_entry_list)

        else:
            intrachrom_entry_end = str(int(annot_sv_obj.start) + 1)
            bedpe_intrachrom_tup = (annot_sv_obj.chrom, annot_sv_obj.start, intrachrom_entry_end,
                                     annot_sv_obj.partner_chrom, annot_sv_obj.partner_start, annot_sv_obj.partner_end,
                                     annot_sv_obj.type, score, annot_sv_obj.strand_1, annot_sv_obj.strand_2,
                                     annot_sv_obj.cytoband, support, annotations)
            high_quality_sv_bedpe.write('{0}\n'.format(bedpe_output_writer(bedpe_intrachrom_tup)))

    # Identify breakend mate records so they can be merged into a single output record
    breakend_with_perfect_mate = {}
    breakend_with_imperfect_mate = []
    for i in range(len(inter_chrom_list)):
        breakend_id_1 = ':'.join((inter_chrom_list[i][0:2]))

        for j in range(len(inter_chrom_list)):
            breakend_id_2 = ':'.join((inter_chrom_list[j][3:5]))

            if breakend_id_1 == breakend_id_2:

                if j in breakend_with_perfect_mate.keys():
                    break
                else:
                    breakend_with_perfect_mate[i] = j

        if i not in breakend_with_perfect_mate.keys() and i not in breakend_with_perfect_mate.values() and i not in breakend_with_imperfect_mate:
            breakend_with_imperfect_mate.append(i)

    # Merge perfect mate breakend record annotations into single output
    for breakend_index in breakend_with_perfect_mate.items():

        base_record = '\t'.join((inter_chrom_list[breakend_index[0]][0:10]))

        end1_cytoband = inter_chrom_list[breakend_index[0]][0].replace('chr','') + inter_chrom_list[breakend_index[0]][10]
        end2_cytoband = inter_chrom_list[breakend_index[1]][0].replace('chr','') + inter_chrom_list[breakend_index[1]][10]
        merged_cytoband = end1_cytoband + ";" + end2_cytoband

        shared_read_support = inter_chrom_list[breakend_index[0]][11]

        end1_original_annotations = inter_chrom_list[breakend_index[0]][12].split(';')
        end2_original_annotations = inter_chrom_list[breakend_index[1]][12].split(';')
        end1_genes = end1_original_annotations[0].replace('=', '_end1=')
        end2_genes = end2_original_annotations[0].replace('=', '_end2=')
        merged_genes = ';'.join((end1_genes, end2_genes))
        end1_gc_content = end1_original_annotations[1]
        end2_gc_content = end2_original_annotations[1].replace('end1', 'end2')
        merged_gc_content = ';'.join((end1_gc_content, end2_gc_content))
        end1_repeat_info = ';'.join(end1_original_annotations[3:5])
        end2_repeat_info = ';'.join(end2_original_annotations[3:5]).replace('end1', 'end2')
        merged_repeat_info = ';'.join((end1_repeat_info, end2_repeat_info))
        end1_gap = end1_original_annotations[7]
        end2_gap = end2_original_annotations[7].replace('end1', 'end2')
        merged_gap = ';'.join((end1_gap, end2_gap))
        end1_segdup = end1_original_annotations[9]
        end2_segdup = end2_original_annotations[9].replace('end1', 'end2')
        merged_segdup = ';'.join((end1_segdup, end2_segdup))
        end1_encode_info = ';'.join(end1_original_annotations[11:13])
        end2_encode_info = ';'.join(end2_original_annotations[11:13]).replace('end1', 'end2')
        merged_encode_info = ';'.join((end1_encode_info, end2_encode_info))
        end1_annotsv_id = end1_original_annotations[15].replace('=', '_end1=')
        end2_annotsv_id = end2_original_annotations[15].replace('=', '_end2=')
        merged_annotsv_id = ';'.join((end1_annotsv_id, end2_annotsv_id))

        merged_annotations = ';'.join((merged_genes, merged_gc_content, merged_repeat_info, merged_gap, merged_segdup, merged_encode_info, merged_annotsv_id))
        merged_interchrom_tup = (base_record, merged_cytoband, shared_read_support, merged_annotations)

        high_quality_sv_bedpe.write('{0}\n'.format(bedpe_output_writer(merged_interchrom_tup)))

    # For breakend records that don't perfectly align with mate end, print both records with annotations of end1
    for imperfect_mate_index in breakend_with_imperfect_mate:

        single_imperfect_breakend_base_record = '\t'.join(inter_chrom_list[imperfect_mate_index][0:10])

        single_imperfect_breakend_cytoband = inter_chrom_list[imperfect_mate_index][0].replace('chr','') + inter_chrom_list[imperfect_mate_index][10]

        single_imperfect_breakend_support = inter_chrom_list[imperfect_mate_index][11]

        single_imperfect_breakend_original_annotations = inter_chrom_list[imperfect_mate_index][12].split(';')
        single_imperfect_breakend_flag = "IMPERFECT"
        single_imperfect_breakend_end1_annotations = ';'.join([single_imperfect_breakend_flag, single_imperfect_breakend_original_annotations[0],
                                                               single_imperfect_breakend_original_annotations[1], single_imperfect_breakend_original_annotations[3],
                                                               single_imperfect_breakend_original_annotations[4], single_imperfect_breakend_original_annotations[7],
                                                               single_imperfect_breakend_original_annotations[9], single_imperfect_breakend_original_annotations[11],
                                                               single_imperfect_breakend_original_annotations[12], single_imperfect_breakend_original_annotations[15]])
        imperfect_breakend_tup = (single_imperfect_breakend_base_record, single_imperfect_breakend_cytoband,
                                  single_imperfect_breakend_support, single_imperfect_breakend_end1_annotations)

        high_quality_sv_bedpe.write('{0}\n'.format(bedpe_output_writer(imperfect_breakend_tup)))

high_quality_sv_bedpe.close()