#!/usr/bin/env Rscript

# This script accepts a fit copy number profile CSV and segmented LogR file from
# Battenberg and generates a CNV profile of chained segments that match the
# clonal nonrounded profile

#########################
#####   Libraries   #####

suppressPackageStartupMessages(library(devgru))
suppressPackageStartupMessages(library(gUtils))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doParallel))

Sys.setenv(DEFAULT_BSGENOME = 'BSgenome.Hsapiens.UCSC.hg38::Hsapiens')
options(scipen = 999)


#########################
#####   Functions   #####

# Battenberg fit.copy.number output currently only provides start point of segment, need to add the end
add_battenberg_segment_end <- function(battenberg_cnv_segments, chrom_iter_list, threads) {
  
  # Set number of threads for foreach loop and file writing
  registerDoParallel(cores = threads)
  
  proper_bed_style_battenberg <- foreach(chrom = 1:length(chrom_iter_list), .combine = rrbind, .packages = "gUtils") %dopar% {
    
    # Subset Battenberg CNV segments by chromosome
    chrom_specific_subset <- battenberg_cnv_segments[seqnames == chrom_iter_list[chrom]]
    
    # Create an empty DT of all segments, now with end column
    complete_segments <- data.table::data.table()
    
    # Start of loop through all records per chrom
    for(i in 1:nrow(chrom_specific_subset)) {
      
      # First, to find the end of each segment we need a record (i) and its successor (i+1)
      segment_slice <- chrom_specific_subset[c(i,i+1)]
      
      # Edge case: for final record of a chromosome, i+1 is all NAs
      # Will account for this case by setting end = start of i-th record
      if(is.na(segment_slice[2, start])) {
        # Edge case
        complete_segments <- gUtils::rrbind(complete_segments,
                                            data.table::data.table(seqnames = segment_slice[1, seqnames],
                                                                   start = segment_slice[1, start],
                                                                   end = segment_slice[1, start],
                                                                   total_cn_nonrounded = segment_slice[1, total_cn_nonrounded],
                                                                   major_allele_nonrounded = segment_slice[1, major_allele_nonrounded],
                                                                   minor_allele_nonrounded = segment_slice[1, minor_allele_nonrounded],
                                                                   segmented_baf = segment_slice[1, segmented_baf],
                                                                   segmented_logr = segment_slice[1, segmented_logr],
                                                                   total_cn_rounded = segment_slice[1, total_cn_rounded],
                                                                   major_allele_rounded = segment_slice[1, major_allele_rounded],
                                                                   minor_allele_rounded = segment_slice[1, minor_allele_rounded]))
        
      } else {
        # For normal records, set end = start of i+1 record
        complete_segments <- gUtils::rrbind(complete_segments,
                                            data.table::data.table(seqnames = segment_slice[1, seqnames],
                                                                   start = segment_slice[1, start],
                                                                   end = segment_slice[2, start],
                                                                   total_cn_nonrounded = segment_slice[1, total_cn_nonrounded],
                                                                   major_allele_nonrounded = segment_slice[1, major_allele_nonrounded],
                                                                   minor_allele_nonrounded = segment_slice[1, minor_allele_nonrounded],
                                                                   segmented_baf = segment_slice[1, segmented_baf],
                                                                   segmented_logr = segment_slice[1, segmented_logr],
                                                                   total_cn_rounded = segment_slice[1, total_cn_rounded],
                                                                   major_allele_rounded = segment_slice[1, major_allele_rounded],
                                                                   minor_allele_rounded = segment_slice[1, minor_allele_rounded]))
      }
    }
    # Print out results of the foreach loop to capture the output
    complete_segments
  }
  # return the final DT
  return(data.table::as.data.table(proper_bed_style_battenberg))
}

# Loop through all CNV segments to chain all segments that have consistent nonrounded major/minor alleles
chain_battenberg_segments <- function(cnv_segments_file, chrom_iter_list, threads) {
  
  # Set number of threads for foreach loop and file writing
  registerDoParallel(cores = threads)
  
  chained_cnv_profile <- foreach(chrom = 1:length(chrom_iter_list), .combine = rrbind, .packages = "gUtils") %dopar% {
    
    # Subset whole CNV dataset by chromosome
    chrom_specific_cnv_segments <- cnv_segments_file[seqnames == chrom_iter_list[chrom]]
    
    # Create a DT of all segments with consistent nonrounded major/minor alleles
    consistent_segments <- data.table::data.table()
    
    # Start of loop through all records per chrom
    for(i in 1:nrow(chrom_specific_cnv_segments)) {
      
      # Add record to start segment chaining if consistent DT is empty, if not check next record
      if(nrow(consistent_segments) == 0) {
        consistent_segments <- gUtils::rrbind(consistent_segments,
                                              data.table::data.table(
                                                seqnames = chrom_specific_cnv_segments[i, seqnames],
                                                start = chrom_specific_cnv_segments[i, start],
                                                end = chrom_specific_cnv_segments[i, end],
                                                total_cn_nonrounded = chrom_specific_cnv_segments[i, total_cn_nonrounded],
                                                major_allele_nonrounded = chrom_specific_cnv_segments[i, major_allele_nonrounded],
                                                minor_allele_nonrounded = chrom_specific_cnv_segments[i, minor_allele_nonrounded],
                                                segmented_baf = chrom_specific_cnv_segments[i, segmented_baf],
                                                segmented_logr = chrom_specific_cnv_segments[i, segmented_logr],
                                                total_cn_rounded = chrom_specific_cnv_segments[i, total_cn_rounded],
                                                major_allele_rounded = chrom_specific_cnv_segments[i, major_allele_rounded],
                                                minor_allele_rounded = chrom_specific_cnv_segments[i, minor_allele_rounded]))
      }
      
      # Check if next record has consistent total copy number and major/minor alleles as previous record
      last_segment_index <- nrow(consistent_segments)
      
      if(chrom_specific_cnv_segments[i, total_cn_nonrounded] == consistent_segments[last_segment_index, "total_cn_nonrounded"] & chrom_specific_cnv_segments[i, major_allele_nonrounded] == consistent_segments[last_segment_index, "major_allele_nonrounded"] & chrom_specific_cnv_segments[i, minor_allele_nonrounded] == consistent_segments[last_segment_index, "minor_allele_nonrounded"]) {
        
        # If segment is consistent with previous, create new record that merges the segment ends
        chained_consistent_segment <- data.table::data.table(
          seqnames = consistent_segments[last_segment_index, "seqnames"],
          start = consistent_segments[last_segment_index, "start"],
          end = chrom_specific_cnv_segments[i, end],
          total_cn_nonrounded = consistent_segments[last_segment_index, "total_cn_nonrounded"],
          major_allele_nonrounded = consistent_segments[last_segment_index, "major_allele_nonrounded"],
          minor_allele_nonrounded = consistent_segments[last_segment_index, "minor_allele_nonrounded"],
          segmented_baf = consistent_segments[last_segment_index, "segmented_baf"],
          segmented_logr = consistent_segments[last_segment_index, "segmented_logr"],
          total_cn_rounded = consistent_segments[last_segment_index, "total_cn_rounded"],
          major_allele_rounded = consistent_segments[last_segment_index, "major_allele_rounded"],
          minor_allele_rounded = consistent_segments[last_segment_index, "minor_allele_rounded"])
        
        # Replace the segment record with new chained segment
        consistent_segments[last_segment_index,] <- chained_consistent_segment
        
      } else {
        # If not consistent, add new segment as a record to then check for chaining
        consistent_segments <- gUtils::rrbind(consistent_segments,
                                              data.table::data.table(
                                                seqnames = chrom_specific_cnv_segments[i, seqnames],
                                                start = chrom_specific_cnv_segments[i, start],
                                                end = chrom_specific_cnv_segments[i, end],
                                                total_cn_nonrounded = chrom_specific_cnv_segments[i, total_cn_nonrounded],
                                                major_allele_nonrounded = chrom_specific_cnv_segments[i, major_allele_nonrounded],
                                                minor_allele_nonrounded = chrom_specific_cnv_segments[i, minor_allele_nonrounded],
                                                segmented_baf = chrom_specific_cnv_segments[i, segmented_baf],
                                                segmented_logr = chrom_specific_cnv_segments[i, segmented_logr],
                                                total_cn_rounded = chrom_specific_cnv_segments[i, total_cn_rounded],
                                                major_allele_rounded = chrom_specific_cnv_segments[i, major_allele_rounded],
                                                minor_allele_rounded = chrom_specific_cnv_segments[i, minor_allele_rounded]))
      }
    }
    # Print out results of the foreach loop to capture the output
    consistent_segments
  }
  # return the final DT
  return(chained_cnv_profile)
}


#########################
#####   Execution   #####

# Accept command line arguments as input
input_args <- commandArgs(trailingOnly = T)

tumor_normal_id <- input_args[1]

fit_cnv_profile_csv <- input_args[2]

segmented_logr_file <- input_args[3]

threads <- input_args[4]

output_bed_filename <- input_args[5]

message("Reading in segmented Battenberg input.....")
cnv_csv <- data.table::fread(input = fit_cnv_profile_csv,
                             sep = ",",
                             header = T,
                             nThread = threads)

seg_logr <- data.table::fread(input = segmented_logr_file,
                              sep = "\t",
                              header = F,
                              col.names = c("seqnames", "start", "segmentedR"),
                              nThread = threads)

cnv_profile_per_segment <- data.table::data.table(seqnames = paste0("chr",seg_logr$seqnames),
                                                  start = seg_logr$start,
                                                  total_cn_nonrounded = cnv_csv$nAfull + cnv_csv$nBfull,
                                                  major_allele_nonrounded = cnv_csv$nAfull,
                                                  minor_allele_nonrounded = cnv_csv$nBfull,
                                                  segmented_baf = cnv_csv$segmentedBAF,
                                                  segmented_logr = cnv_csv$segmentedR,
                                                  total_cn_rounded = cnv_csv$nA + cnv_csv$nB,
                                                  major_allele_rounded = cnv_csv$nA,
                                                  minor_allele_rounded = cnv_csv$nB)

# Split each input loop by chromosome
chrom_loop_list <- as.character(unique(cnv_profile_per_segment$seqnames))

message("Adding proper end to each CNV segment.....")
complete_end_to_end_segments <- add_battenberg_segment_end(battenberg_cnv_segments = cnv_profile_per_segment,
                                                           chrom_iter_list = chrom_loop_list,
                                                           threads = threads)

message("Chaining CNV segments that have consistent nonrounded total CN and major/minor alleles.....")
chained_cnv_alleles <- chain_battenberg_segments(cnv_segments_file = complete_end_to_end_segments,
                                                 chrom_iter_list = chrom_loop_list,
                                                 threads = threads)

message("Writing output CNV profile BED.....")
data.table::fwrite(x = chained_cnv_alleles,
                   file = ,
                   sep = "\t",
                   row.names = F,
                   col.names = T)
message("D O N E .....")

