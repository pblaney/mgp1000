#!/usr/bin/env Rscript

# This script accepts a text file containing per variant REF and ALT supporting
# reads on the forward and reverse strands to perform a simple test of minimum
# proportion of reads in both directions.
# It then produces a BED file with each variant's coordinates that fail the test
# to be used in filtering.

#########################
#####   Libraries   #####

suppressPackageStartupMessages(library(tidyverse))

#########################
#####   Functions   #####

strand_bias_proportion_tester <- function(strand_metrics_file) {
  
  # Read in input file containing REF/ALT read data per forward and reverse strand
  print("Reading input strand metrics file.....")
  strand_metrics <- read_delim(file = strand_metrics_file,
                               delim = "\t",
                               col_names = c("chrom", "start", "forward_strand", "reverse_strand"),
                               col_types = "cicc")
  
  # Create output data frame that holds BED lines
  strand_bias_filter_bed_file <- tibble()
  
  # Loop through each entry in the strand metrics file and perform test of minimum proportion
  print("Performing test of minimum proportion per entry.....")
  for(i in 1:dim(strand_metrics)[1]) {
    
    # Assign each set of strand reads
    forward_strand_reads <- as.numeric(str_split(string = strand_metrics$forward_strand[i],
                                                 pattern = ",",
                                                 simplify = TRUE))
    reverse_strand_reads <- as.numeric(str_split(string = strand_metrics$reverse_strand[i],
                                                 pattern = ",",
                                                 simplify = TRUE))
    
    # Calculate the minimum proportion of reads in each direction
    strand_read_tibble <- tibble("forward_strand_reads" = forward_strand_reads,
                                 "reverse_strand_reads" = reverse_strand_reads)
    forward_proportion <- sum(strand_read_tibble$forward_strand_reads) / sum(strand_read_tibble)
    reverse_proportion <- sum(strand_read_tibble$reverse_strand_reads) / sum(strand_read_tibble)
    minimum_sb_proportion <- min(forward_proportion, reverse_proportion)
    
    # Populate tibble with evaluation of each strand bias proportion test if the min strand bias proportion <0.10
    if(minimum_sb_proportion < 0.10) {

      strand_bias_filter_bed_file <- rbind(strand_bias_filter_bed_file,
                                           tibble("chrom" = strand_metrics$chrom[i],
                                                  "start" = strand_metrics$start[i],
                                                  "end" = strand_metrics$start[i],
                                                  "minimum_sb_proportion" = minimum_sb_proportion,
                                                  "forward_reads_ref" = forward_strand_reads[1],
                                                  "forward_reads_alt" = forward_strand_reads[2],
                                                  "reverse_reads_ref" = reverse_strand_reads[1],
                                                  "reverse_reads_alt" = reverse_strand_reads[2]))
    }
  }
  
  # Write out BED line for variant to be filtered out based on minimum strand bias proportion
  print("Printing output BED file for filtering.....")
  return(strand_bias_filter_bed_file)
}

#########################
#####   Execution   #####

# Accept command line arguments as input
input_args <- commandArgs(trailingOnly = T)

# Execute function
output_bed <- strand_bias_proportion_tester(input_args[1])
output_bed_filename <- input_args[2]

write.table(x = output_bed,
            file = output_bed_filename,
            sep = "\t",
            col.names = TRUE,
            quote = FALSE,
            row.names = FALSE)
