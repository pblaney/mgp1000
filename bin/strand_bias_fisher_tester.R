#!/usr/bin/env Rscript

# This script accepts a text file containing per variant REF and ALT supporting
# reads on the forward and reverse strands to perform a Fisher's Exact Test.
# It then produces a BED file with each variant's coordinates that fail the test
# to be used in filtering.

suppressPackageStartupMessages(library(tidyverse))

strand_bias_fisher_tester <- function(strand_metrics_file, output_bed_filename) {
  
  # Read in input file containing REF/ALT read data per forward and reverse strand
  print("Reading input strand metrics file.....")
  strand_metrics <- read_delim(file = strand_metrics_file,
                               delim = "\t",
                               col_names = c("chrom", "start", "forward_strand", "reverse_strand"))
  
  # Create output data frame that holds BED lines
  strand_bias_filter_bed_file <- tibble("chrom" = rep(NA, dim(strand_metrics)[1]),
                                        "start" = rep(NA, dim(strand_metrics)[1]),
                                        "end" = rep(NA, dim(strand_metrics)[1]),
                                        "p_value" = rep(NA, dim(strand_metrics)[1]),
                                        "forward_reads_ref" = rep(NA, dim(strand_metrics)[1]),
                                        "forward_reads_alt" = rep(NA, dim(strand_metrics)[1]),
                                        "reverse_reads_ref" = rep(NA, dim(strand_metrics)[1]),
                                        "reverse_reads_alt" = rep(NA, dim(strand_metrics)[1]))
  
  # Loop through each entry in the strand metrics file and perform Fisher's Exact Test
  print("Performing Fisher's Exact Test per entry.....")
  for(i in 1:dim(strand_metrics)[1]) {
    
    # Assign each set of strand reads
    forward_strand_reads <- as.numeric(str_split(string = strand_metrics$forward_strand[i],
                                                 pattern = ",",
                                                 simplify = TRUE))
    reverse_strand_reads <- as.numeric(str_split(string = strand_metrics$reverse_strand[i],
                                                 pattern = ",",
                                                 simplify = TRUE))
    
    # Set up 2x2 contingency table for Fisher's Exact Test
    strand_read_tibble <- tibble("forward_strand_reads" = forward_strand_reads,
                                 "reverse_strand_reads" = reverse_strand_reads)
    rownames(strand_read_tibble) <- c("REF", "ALT")
    
    # Extract p-value
    p_value <- fisher.test(x = as.matrix(strand_read_tibble),
                           alternative = "two.sided")$p.value
    
    # Populate tibble with evaluation of each Fisher's Exact Test
    strand_bias_filter_bed_file$chrom[i] <- strand_metrics$chrom[i]
    strand_bias_filter_bed_file$start[i] <- strand_metrics$start[i]
    strand_bias_filter_bed_file$end[i] <- strand_metrics$start[i]
    strand_bias_filter_bed_file$p_value[i] <- p_value
    strand_bias_filter_bed_file$forward_reads_ref[i] <- forward_strand_reads[1]
    strand_bias_filter_bed_file$forward_reads_alt[i] <- forward_strand_reads[2]
    strand_bias_filter_bed_file$reverse_reads_ref[i] <- reverse_strand_reads[1]
    strand_bias_filter_bed_file$reverse_reads_alt[i] <- reverse_strand_reads[2]
    
  }
  
  # Write out BED line for variant to be filtered out based on  p-value < 0.10
  print("Printing output BED file for filtering.....")
  strand_bias_filter_bed_file %>%
    filter(p_value < 0.1) %>%
    write_delim(file = output_bed_filename,
                delim = "\t",
                col_names = T)
}

# Accept command line arguments as input
input_args <- commandArgs(trailingOnly = T)

# Execute function
strand_bias_fisher_tester(input_args[1], input_args[2])
