#!/usr/bin/env Rscript

# This script accepts consensus SV BEDPE and IgCaller oncogenic rearrangement TSV files
# as input and generates a merged, deduplicated, and simplified BEDPE for easy
# downstream analysis

#########################
#####   Libraries   #####

suppressPackageStartupMessages(library(gGnome))
suppressPackageStartupMessages(library(gUtils))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(tidyverse))

options(scipen = 999)

#########################
#####   Functions   #####

# Modification of the gUtils function gr.collapse, 
gr_collapse = function(gr, pad, strand_ignore) {
  tmp = gr.findoverlaps(gr + pad, gr + pad, ignore.strand = strand_ignore)
  m = rep(FALSE, length(gr))
  m[tmp$query.id[tmp$query.id == (tmp$subject.id-1)]] = TRUE
  
  ## will not collapse if two intersecting ranges are in the wrong "order" (ie not increasing (decreasing) on pos (neg) strand
  m[which((strand(gr)[-length(gr)] == '+' & (start(gr)[-length(gr)] > start(gr)[-1])) |
            (strand(gr)[-length(gr)] == '-' & (end(gr)[-length(gr)] < end(gr)[-1])))] = FALSE
  
  m = as(m, 'IRanges')
  
  if (length(m)>0)
  {
    end(m) = end(m)+1
    tmp = cbind(start(gr)[start(m)], end(gr)[start(m)], start(gr)[end(m)], end(gr)[end(m)])
    s = pmin(start(gr)[start(m)], end(gr)[start(m)], start(gr)[end(m)], end(gr)[end(m)])
    e = pmax(start(gr)[start(m)], end(gr)[start(m)], start(gr)[end(m)], end(gr)[end(m)])
    return(GRanges(seqnames(gr)[start(m)], IRanges(s, e), strand = strand(gr)[start(m)], seqlengths = seqlengths(gr)))
  }
  else{
    return(gr[c()])
  }
}

# Main polishing function
merge_dedup_simplify_svs <- function(tum_norm_id, bedpe, igcaller_onco_tsv) {
  
  # Read in IgCaller oncogeneic rearrangements TSV, filter then convert to BEDPE
  igcaller <- read_delim(file = igcaller_onco_tsv,
                         delim = "\t",
                         show_col_types = FALSE,
                         col_names = T)
  
  igcaller_bedpe <- igcaller %>%
    dplyr::filter(Score >= 5 & `Reads in normal` <= 3 & `Count in PoN` <= 4) %>%
    dplyr::mutate("end1" = PositionA + 1,
                  "end2" = PositionB + 1,
                  "score" = ".") %>%
    dplyr::select(ChrA, PositionA, end1, ChrB, PositionB, end2, Mechanism, score, StrandA, StrandB) %>%
    dplyr::rename("#chr1" = ChrA,
                  "start1" = PositionA,
                  "chr2" = ChrB,
                  "start2" = PositionB,
                  "type" = Mechanism,
                  "strand1" = StrandA,
                  "strand2" = StrandB) %>%
    dplyr::mutate("caller" = "igcaller")
  
  igcaller_bedpe$type <- case_when(igcaller_bedpe$type == "Deletion" ~ "DEL",
                                   igcaller_bedpe$type == "Duplication" ~ "DUP",
                                   igcaller_bedpe$type == "Insertion" ~ "INS",
                                   igcaller_bedpe$type == "Inversion" ~ "INV",
                                   igcaller_bedpe$type == "Translocation" ~ "BND")
  
  # Read in MGP1000 Consensus SV BEDPE
  consensus_bedpe <- read_delim(file = bedpe,
                                delim = "\t",
                                col_names = T,
                                show_col_types = FALSE,
                                col_select = c(1:10, "score")) %>%
    dplyr::mutate("caller" = str_remove(string = as.character(str_split(string = score, pattern = ";", simplify = T)[,1]), 
                                        pattern = "Callers="))
  consensus_bedpe$score <- "."

  # Rename first column if in '#chrom' format
  if(sum(str_detect(string = colnames(consensus_bedpe), pattern = "chrom")) > 0) {

    colnames(consensus_bedpe) <- str_replace(string = colnames(consensus_bedpe), pattern = "chrom", replacement = "chr")
  }

  # Write the simple BEDPEs as temporary files
  temp_igcaller_bedpe <- tempfile(pattern = tum_norm_id, fileext = ".igcaller.simple.bedpe")
  write.table(x = igcaller_bedpe, file = temp_igcaller_bedpe,
              row.names = F, col.names = T, sep = "\t", quote = F)

  temp_consensus_bedpe <- tempfile(pattern = tum_norm_id, fileext = ".consensus.simple.bedpe")
  write.table(x = consensus_bedpe, file = temp_consensus_bedpe,
              row.names = F, col.names = T, sep = "\t", quote = F)
  
  # Read in BEDPE files into gGnome junction objects
  ref_seq_info <- c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
                    159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
                    114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
                    58617616, 64444167, 46709983, 50818468, 156040895, 57227415)
  names(ref_seq_info) <- paste0("chr", c(seq(1,22,1), "X", "Y"))
  
  # Edge case: No records in IgCaller pass filter so this BEDPE is empty
  if(nrow(igcaller_bedpe) == 0) {
    igcaller_jnc <- GRanges()
  } else {
    igcaller_jnc <- jJ(rafile = temp_igcaller_bedpe,
                       chr.convert = F,
                       keep.features = T,
                       seqlengths = ref_seq_info,
                       seqlevels = ref_seq_info)
  }
  
  consensus_jnc <- jJ(rafile = temp_consensus_bedpe,
                      chr.convert = F,
                      keep.features = T,
                      seqlengths = ref_seq_info,
                      seqlevels = ref_seq_info)
  
  merged_jnc <- merge.Junction(igcaller = igcaller_jnc,
                               consensus = consensus_jnc)
  
  # Convert GRangeList of junctions to GRanges
  merged_jnc_gr <- gUtils::grl.unlist(grl = merged_jnc$grl)
  
  # Loop through chromosomes in merged junction GRanges and update the seqlengths
  for(i in 1:length(seqlevels(merged_jnc_gr))) {
    seqlengths(merged_jnc_gr)[seqlevels(merged_jnc_gr)[i]] <- ref_seq_info[seqlevels(merged_jnc_gr)[i]]
  }
  
  # Reset order of chromosomes to standard
  merged_jnc_gr@seqnames@values <- factor(merged_jnc_gr@seqnames@values,
                                          levels = mixedsort(seqlevels(seqinfo(merged_jnc_gr))))
  merged_jnc_gr@seqinfo <- sortSeqlevels(merged_jnc_gr@seqinfo, X.is.sexchrom = T)
  
  merged_jnc_gr <- sort.GenomicRanges(merged_jnc_gr, ignore.strand = TRUE)
  
  # Edge case (continued): code expects igcaller columns based on the merging, create dummy variables if igcaller empty
  if(nrow(igcaller_bedpe) == 0) {
    merged_jnc_gr$name.igcaller <- NA
    merged_jnc_gr$NA..igcaller <- NA
    
    merged_jnc_gr$name.consensus <- merged_jnc_gr$name
    merged_jnc_gr$NA..consensus <- merged_jnc_gr$NA.
  }
  
  # Find all duplicated inter- and intra- chromosomal SVs and set them as ranges to remove
  inter_to_remove <- gr_collapse(merged_jnc_gr[which(merged_jnc_gr$name.consensus == "BND" | merged_jnc_gr$name.igcaller == "BND"),],
                                 pad = 10, strand_ignore = T)
  intra_to_remove <- gr_collapse(merged_jnc_gr[which(merged_jnc_gr$name.consensus != "BND" | merged_jnc_gr$name.igcaller != "BND"),],
                                 pad = 10, strand_ignore = T)
  ranges_to_remove <- grbind(inter_to_remove, intra_to_remove)
  #ranges_to_remove <- gr_collapse(merged_jnc_gr, pad = 10, strand_ignore = T)
  
  # Loop through all SV pairs and find duplicates
  non_duplicate_ranges <- GenomicRanges::GRanges()
  
  # (good) Edge case: No duplicated ranges, just add row to output
  if(is.null(ranges_to_remove)) {
    non_duplicate_ranges <- merged_jnc_gr
    
  } else {
    for (i in 1:n_distinct(merged_jnc_gr$merged.ix)) {
      
      # Get breakpoints of SV pair
      sv_pair <- merged_jnc_gr %Q% (merged.ix == i)
      
      # First, check if either breakpoint in the SV pair is found in the set of ranges to remove,
      # If not, add the good record to the non_duplicated range set
      if(sum(sv_pair %^% ranges_to_remove) == 0) {
        
        non_duplicate_ranges <- grbind(x = non_duplicate_ranges, sv_pair)
        
        # Second, if either of the breakpoints are found in the set of ranges to remove
        # check if the record has already been added to the final output
        # If not, add the record, otherwise continue to next SV pair
      } else if(sum(sv_pair %^% ranges_to_remove) > 0 & sum((sv_pair + 10) %^% non_duplicate_ranges) == 0) {
        
        non_duplicate_ranges <- grbind(x = non_duplicate_ranges, sv_pair)
      }
    }
  }
  
  # Deduplicated SV GRanges
  dedup_jnc_gr <- merged_jnc_gr[which(merged_jnc_gr$merged.ix %in% non_duplicate_ranges$merged.ix),]
  
  # Convert GRanges to dataframe for easy manipulation
  dedup_jnc_df <- gr2dt(dedup_jnc_gr)
  
  # Construct the final BEDPE output
  final_bedpe <- tibble()
  for (i in 1:length(unique(dedup_jnc_df$merged.ix))) {
    
    # Get breakpoints of SV pair
    bedpe_pair <- dedup_jnc_df %>%
      dplyr::filter(merged.ix == unique(dedup_jnc_df$merged.ix)[i]) %>%
      dplyr::select(seqnames, start, end, strand, name.igcaller, name.consensus, NA..igcaller, NA..consensus)
    
    record_type <- na.omit(unique(c(bedpe_pair$name.igcaller, bedpe_pair$name.consensus)))[1]
    
    record_caller <- na.omit(unique(c(bedpe_pair$NA..igcaller, bedpe_pair$NA..consensus)))[1]
    
    bedpe_single_line_record <- tibble("#chr1" = bedpe_pair$seqnames[1],
                                       "start1" = bedpe_pair$start[1],
                                       "end1" = bedpe_pair$end[1],
                                       "chr2" = bedpe_pair$seqnames[2],
                                       "start2" = bedpe_pair$start[2],
                                       "end2" = bedpe_pair$end[2],
                                       "type" = record_type,
                                       "score" = ".",
                                       "strand1" = bedpe_pair$strand[1],
                                       "strand2" = bedpe_pair$strand[2],
                                       "caller" = record_caller)
    
    final_bedpe <- rbind(final_bedpe,
                         bedpe_single_line_record)
  }
  
  # Remove temp files and return output
  unlink(temp_igcaller_bedpe)
  unlink(temp_consensus_bedpe)
  return(final_bedpe)
}

#########################
#####   Execution   #####

# Accept command line arguments as input
input_args <- commandArgs(trailingOnly = T)

sv_file_dir <- input_args[1]

output_dir <- input_args[2]

# Find all input files in directory provided
print("Scanning directory for MGP1000 consensus SV BEDPE and IgCaller oncogenic rearrangment TSV files.....")
igcaller_tsv_input <- list.files(path = sv_file_dir,
                                 pattern = "*.igcaller.oncogenic.rearrangements.tsv")

consensus_bedpe_input <- list.files(path = sv_file_dir,
                                    pattern = "*.hq.consensus.somatic.sv.annotated.bedpe")

# Get sample names and gather input file pairs
print("Collecting possible input file pairs.....")
igcaller_samples <- str_remove(string = igcaller_tsv_input, pattern = ".igcaller.oncogenic.rearrangements.tsv") %>%
                      as_tibble_col(column_name = "sample")

consensus_samples <- str_remove(string = consensus_bedpe_input, pattern = ".hq.consensus.somatic.sv.annotated.bedpe") %>%
                       as_tibble_col(column_name = "sample")

sample_set <- inner_join(x = igcaller_samples, y = consensus_samples)

# Generate final BEDPE from input pair
print("Generating final merged, deduplicated, simplified SV BEDPE for possible input pairs.....")
for(i in 1:nrow(sample_set)) {
  igcaller_tsv <- igcaller_tsv_input[str_detect(string = igcaller_tsv_input, pattern = sample_set$sample[i]) %>% which()]
  
  consensus_bedpe <- consensus_bedpe_input[str_detect(string = consensus_bedpe_input, pattern = sample_set$sample[i]) %>% which()]
  
  # Run main polishing function 
  final_bedpe_df <- merge_dedup_simplify_svs(tum_norm_id = sample_set$sample[i],
                                             bedpe = paste0(sv_file_dir, consensus_bedpe),
                                             igcaller_onco_tsv = paste0(sv_file_dir, igcaller_tsv))
  
  # Write final BEDPE output file
  polished_bedpe_filename <- paste0(sample_set$sample[i], ".hq.final.consensus.somatic.sv.bedpe")
  
  write.table(x = final_bedpe_df,
              file = paste0(output_dir, polished_bedpe_filename),
              sep = "\t",
              col.names = TRUE,
              quote = FALSE,
              row.names = FALSE)

  print(paste(polished_bedpe_filename, "..... D O N E"))
}
