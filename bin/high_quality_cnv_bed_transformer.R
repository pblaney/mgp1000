#!/usr/bin/env Rscript

# This script accepts a BED file of merged consensus CNV and alleles per segment
# and a string to set the sex of the sample.
# It then produces a concise high quality BED file with chained segments and
# robust classification

#########################
#####   Libraries   #####

suppressPackageStartupMessages(library(tidyverse))

#########################
#####   Functions   #####

segment_disagreement_resolution <- function(consensus_cnv_alleles_bed_df, consensus_mechanism) {

  agreement_resolved_cnv_df <-  tibble()

  # First, check which consensus mechanism was used: three_way or four_way
  if(consensus_mechanism == "four_way") {

    # Resolve all no agreement segments using Control-FREEC, Sclust, Accucopy, then ascatNgs hierarchy
    for(i in 1:nrow(consensus_cnv_alleles_bed_df)) {
  
      # First find segments with full disagreement
      if(consensus_cnv_alleles_bed_df[i,]$caller_agreement == "no_agreement" &
         consensus_cnv_alleles_bed_df[i,]$allele_caller_agreement == "no_agreement") {
    
        if(consensus_cnv_alleles_bed_df[i,]$controlfreec_cn != ".") {

          control_freec_major_allele <- str_split(string = consensus_cnv_alleles_bed_df[i,]$controlfreec_alleles,
                                                  pattern = "/",
                                                  simplify = T)[1]
          control_freec_minor_allele <- str_split(string = consensus_cnv_alleles_bed_df[i,]$controlfreec_alleles,
                                                  pattern = "/",
                                                  simplify = T)[2]

          agreement_resolved_cnv_df <- rbind(agreement_resolved_cnv_df,
                                             tibble("chrom" = consensus_cnv_alleles_bed_df[i,]$chrom,
                                                    "start" = consensus_cnv_alleles_bed_df[i,]$start,
                                                    "end" = consensus_cnv_alleles_bed_df[i,]$end,
                                                    "consensus_total_cn" = consensus_cnv_alleles_bed_df[i,]$controlfreec_cn,
                                                    "consensus_major_allele" = control_freec_major_allele,
                                                    "consensus_minor_allele" = control_freec_minor_allele,
                                                    "caller_agreement" = "controlfreec",
                                                    "allele_caller_agreement" = "controlfreec"))

        } else if(consensus_cnv_alleles_bed_df[i,]$sclust_cn != ".") {

          sclust_major_allele <- str_split(string = consensus_cnv_alleles_bed_df[i,]$sclust_alleles,
                                           pattern = "/",
                                           simplify = T)[1]
          sclust_minor_allele <- str_split(string = consensus_cnv_alleles_bed_df[i,]$sclust_alleles,
                                           pattern = "/",
                                           simplify = T)[2]
          agreement_resolved_cnv_df <- rbind(agreement_resolved_cnv_df,
                                             tibble("chrom" = consensus_cnv_alleles_bed_df[i,]$chrom,
                                                    "start" = consensus_cnv_alleles_bed_df[i,]$start,
                                                    "end" = consensus_cnv_alleles_bed_df[i,]$end,
                                                    "consensus_total_cn" = consensus_cnv_alleles_bed_df[i,]$sclust_cn,
                                                    "consensus_major_allele" = sclust_major_allele,
                                                    "consensus_minor_allele" = sclust_minor_allele,
                                                    "caller_agreement" = "sclust",
                                                    "allele_caller_agreement" = "sclust"))

        } else if(consensus_cnv_alleles_bed_df[i,]$accucopy_cn != ".") {

          accucopy_major_allele <- str_split(string = consensus_cnv_alleles_bed_df[i,]$accucopy_alleles,
                                             pattern = "/",
                                             simplify = T)[1]
          accucopy_minor_allele <- str_split(string = consensus_cnv_alleles_bed_df[i,]$accucopy_alleles,
                                             pattern = "/",
                                             simplify = T)[2]
          agreement_resolved_cnv_df <- rbind(agreement_resolved_cnv_df,
                                             tibble("chrom" = consensus_cnv_alleles_bed_df[i,]$chrom,
                                                    "start" = consensus_cnv_alleles_bed_df[i,]$start,
                                                    "end" = consensus_cnv_alleles_bed_df[i,]$end,
                                                    "consensus_total_cn" = consensus_cnv_alleles_bed_df[i,]$accucopy_cn,
                                                    "consensus_major_allele" = accucopy_major_allele,
                                                    "consensus_minor_allele" = accucopy_minor_allele,
                                                    "caller_agreement" = "accucopy",
                                                    "allele_caller_agreement" = "accucopy"))
        }
    
        # Then find all segments that only disagree on allele balance, likely due to called LOH by Control-FREEC  
      } else if(consensus_cnv_alleles_bed_df[i,]$caller_agreement != "no_agreement" &
                consensus_cnv_alleles_bed_df[i,]$allele_caller_agreement == "no_agreement") {
    
        if(str_detect(string = consensus_cnv_alleles_bed_df[i,]$caller_agreement, pattern = "controlfreec") &
           consensus_cnv_alleles_bed_df[i,]$controlfreec_cn != ".") {

          control_freec_major_allele <- str_split(string = consensus_cnv_alleles_bed_df[i,]$controlfreec_alleles,
                                                  pattern = "/",
                                                  simplify = T)[1]
          control_freec_minor_allele <- str_split(string = consensus_cnv_alleles_bed_df[i,]$controlfreec_alleles,
                                                  pattern = "/",
                                                  simplify = T)[2]
      
          agreement_resolved_cnv_df <- rbind(agreement_resolved_cnv_df,
                                             tibble("chrom" = consensus_cnv_alleles_bed_df[i,]$chrom,
                                                    "start" = consensus_cnv_alleles_bed_df[i,]$start,
                                                    "end" = consensus_cnv_alleles_bed_df[i,]$end,
                                                    "consensus_total_cn" = consensus_cnv_alleles_bed_df[i,]$consensus_cn,
                                                    "consensus_major_allele" = control_freec_major_allele,
                                                    "consensus_minor_allele" = control_freec_minor_allele,
                                                    "caller_agreement" = consensus_cnv_alleles_bed_df[i,]$caller_agreement,
                                                    "allele_caller_agreement" = "controlfreec"))
      
        } else if(str_detect(string = consensus_cnv_alleles_bed_df[i,]$caller_agreement, pattern = "sclust") &
                  consensus_cnv_alleles_bed_df[i,]$sclust_cn != ".") {

          sclust_major_allele <- str_split(string = consensus_cnv_alleles_bed_df[i,]$sclust_alleles,
                                           pattern = "/",
                                           simplify = T)[1]
          sclust_minor_allele <- str_split(string = consensus_cnv_alleles_bed_df[i,]$sclust_alleles,
                                           pattern = "/",
                                           simplify = T)[2]
          agreement_resolved_cnv_df <- rbind(agreement_resolved_cnv_df,
                                             tibble("chrom" = consensus_cnv_alleles_bed_df[i,]$chrom,
                                                    "start" = consensus_cnv_alleles_bed_df[i,]$start,
                                                    "end" = consensus_cnv_alleles_bed_df[i,]$end,
                                                    "consensus_total_cn" = consensus_cnv_alleles_bed_df[i,]$consensus_cn,
                                                    "consensus_major_allele" = sclust_major_allele,
                                                    "consensus_minor_allele" = sclust_minor_allele,
                                                    "caller_agreement" = consensus_cnv_alleles_bed_df[i,]$caller_agreement,
                                                    "allele_caller_agreement" = "sclust"))

        } else if(str_detect(string = consensus_cnv_alleles_bed_df[i,]$caller_agreement, pattern = "accucopy") &
                  consensus_cnv_alleles_bed_df[i,]$accucopy_cn != ".") {

          accucopy_major_allele <- str_split(string = consensus_cnv_alleles_bed_df[i,]$accucopy_alleles,
                                             pattern = "/",
                                             simplify = T)[1]
          accucopy_minor_allele <- str_split(string = consensus_cnv_alleles_bed_df[i,]$accucopy_alleles,
                                             pattern = "/",
                                             simplify = T)[2]
          agreement_resolved_cnv_df <- rbind(agreement_resolved_cnv_df,
                                             tibble("chrom" = consensus_cnv_alleles_bed_df[i,]$chrom,
                                                    "start" = consensus_cnv_alleles_bed_df[i,]$start,
                                                    "end" = consensus_cnv_alleles_bed_df[i,]$end,
                                                    "consensus_total_cn" = consensus_cnv_alleles_bed_df[i,]$consensus_cn,
                                                    "consensus_major_allele" = accucopy_major_allele,
                                                    "consensus_minor_allele" = accucopy_minor_allele,
                                                    "caller_agreement" = consensus_cnv_alleles_bed_df[i,]$caller_agreement,
                                                    "allele_caller_agreement" = "accucopy"))
        }

      } else {

        agreement_resolved_cnv_df <- rbind(agreement_resolved_cnv_df,
                                           tibble("chrom" = consensus_cnv_alleles_bed_df[i,]$chrom,
                                                  "start" = consensus_cnv_alleles_bed_df[i,]$start,
                                                  "end" = consensus_cnv_alleles_bed_df[i,]$end,
                                                  "consensus_total_cn" = consensus_cnv_alleles_bed_df[i,]$consensus_cn,
                                                  "consensus_major_allele" = consensus_cnv_alleles_bed_df[i,]$consensus_major_allele,
                                                  "consensus_minor_allele" = consensus_cnv_alleles_bed_df[i,]$consensus_minor_allele,
                                                  "caller_agreement" = consensus_cnv_alleles_bed_df[i,]$caller_agreement,
                                                  "allele_caller_agreement" = consensus_cnv_alleles_bed_df[i,]$allele_caller_agreement))
      }
    }

  } else if(consensus_mechanism == "three_way") {

    # Resolve all no agreement segments using Control-FREEC, Accucopy, then ascatNgs hierarchy
    for(i in 1:nrow(consensus_cnv_alleles_bed_df)) {

      # First find segments with full disagreement
      if(consensus_cnv_alleles_bed_df[i,]$caller_agreement == "no_agreement" &
         consensus_cnv_alleles_bed_df[i,]$allele_caller_agreement == "no_agreement") {
    
        if(consensus_cnv_alleles_bed_df[i,]$controlfreec_cn != ".") {

          control_freec_major_allele <- str_split(string = consensus_cnv_alleles_bed_df[i,]$controlfreec_alleles,
                                                  pattern = "/",
                                                  simplify = T)[1]
          control_freec_minor_allele <- str_split(string = consensus_cnv_alleles_bed_df[i,]$controlfreec_alleles,
                                                  pattern = "/",
                                                  simplify = T)[2]

          agreement_resolved_cnv_df <- rbind(agreement_resolved_cnv_df,
                                             tibble("chrom" = consensus_cnv_alleles_bed_df[i,]$chrom,
                                                    "start" = consensus_cnv_alleles_bed_df[i,]$start,
                                                    "end" = consensus_cnv_alleles_bed_df[i,]$end,
                                                    "consensus_total_cn" = consensus_cnv_alleles_bed_df[i,]$controlfreec_cn,
                                                    "consensus_major_allele" = control_freec_major_allele,
                                                    "consensus_minor_allele" = control_freec_minor_allele,
                                                    "caller_agreement" = "controlfreec",
                                                    "allele_caller_agreement" = "controlfreec"))

        } else if(consensus_cnv_alleles_bed_df[i,]$accucopy_cn != ".") {

          accucopy_major_allele <- str_split(string = consensus_cnv_alleles_bed_df[i,]$accucopy_alleles,
                                             pattern = "/",
                                             simplify = T)[1]
          accucopy_minor_allele <- str_split(string = consensus_cnv_alleles_bed_df[i,]$accucopy_alleles,
                                             pattern = "/",
                                             simplify = T)[2]
          agreement_resolved_cnv_df <- rbind(agreement_resolved_cnv_df,
                                             tibble("chrom" = consensus_cnv_alleles_bed_df[i,]$chrom,
                                                    "start" = consensus_cnv_alleles_bed_df[i,]$start,
                                                    "end" = consensus_cnv_alleles_bed_df[i,]$end,
                                                    "consensus_total_cn" = consensus_cnv_alleles_bed_df[i,]$accucopy_cn,
                                                    "consensus_major_allele" = accucopy_major_allele,
                                                    "consensus_minor_allele" = accucopy_minor_allele,
                                                    "caller_agreement" = "accucopy",
                                                    "allele_caller_agreement" = "accucopy"))
        }
    
        # Then find all segments that only disagree on allele balance, likely due to called LOH by Control-FREEC  
      } else if(consensus_cnv_alleles_bed_df[i,]$caller_agreement != "no_agreement" &
                consensus_cnv_alleles_bed_df[i,]$allele_caller_agreement == "no_agreement") {
    
        if(str_detect(string = consensus_cnv_alleles_bed_df[i,]$caller_agreement, pattern = "controlfreec") &
           consensus_cnv_alleles_bed_df[i,]$controlfreec_cn != ".") {

          control_freec_major_allele <- str_split(string = consensus_cnv_alleles_bed_df[i,]$controlfreec_alleles,
                                                  pattern = "/",
                                                  simplify = T)[1]
          control_freec_minor_allele <- str_split(string = consensus_cnv_alleles_bed_df[i,]$controlfreec_alleles,
                                                  pattern = "/",
                                                  simplify = T)[2]
      
          agreement_resolved_cnv_df <- rbind(agreement_resolved_cnv_df,
                                             tibble("chrom" = consensus_cnv_alleles_bed_df[i,]$chrom,
                                                    "start" = consensus_cnv_alleles_bed_df[i,]$start,
                                                    "end" = consensus_cnv_alleles_bed_df[i,]$end,
                                                    "consensus_total_cn" = consensus_cnv_alleles_bed_df[i,]$consensus_cn,
                                                    "consensus_major_allele" = control_freec_major_allele,
                                                    "consensus_minor_allele" = control_freec_minor_allele,
                                                    "caller_agreement" = consensus_cnv_alleles_bed_df[i,]$caller_agreement,
                                                    "allele_caller_agreement" = "controlfreec"))
      
        } else if(str_detect(string = consensus_cnv_alleles_bed_df[i,]$caller_agreement, pattern = "accucopy") &
                  consensus_cnv_alleles_bed_df[i,]$accucopy_cn != ".") {

          accucopy_major_allele <- str_split(string = consensus_cnv_alleles_bed_df[i,]$accucopy_alleles,
                                             pattern = "/",
                                             simplify = T)[1]
          accucopy_minor_allele <- str_split(string = consensus_cnv_alleles_bed_df[i,]$accucopy_alleles,
                                             pattern = "/",
                                             simplify = T)[2]
          agreement_resolved_cnv_df <- rbind(agreement_resolved_cnv_df,
                                             tibble("chrom" = consensus_cnv_alleles_bed_df[i,]$chrom,
                                                    "start" = consensus_cnv_alleles_bed_df[i,]$start,
                                                    "end" = consensus_cnv_alleles_bed_df[i,]$end,
                                                    "consensus_total_cn" = consensus_cnv_alleles_bed_df[i,]$consensus_cn,
                                                    "consensus_major_allele" = accucopy_major_allele,
                                                    "consensus_minor_allele" = accucopy_minor_allele,
                                                    "caller_agreement" = consensus_cnv_alleles_bed_df[i,]$caller_agreement,
                                                    "allele_caller_agreement" = "accucopy"))
        }
      } else {

        agreement_resolved_cnv_df <- rbind(agreement_resolved_cnv_df,
                                           tibble("chrom" = consensus_cnv_alleles_bed_df[i,]$chrom,
                                                  "start" = consensus_cnv_alleles_bed_df[i,]$start,
                                                  "end" = consensus_cnv_alleles_bed_df[i,]$end,
                                                  "consensus_total_cn" = consensus_cnv_alleles_bed_df[i,]$consensus_cn,
                                                  "consensus_major_allele" = consensus_cnv_alleles_bed_df[i,]$consensus_major_allele,
                                                  "consensus_minor_allele" = consensus_cnv_alleles_bed_df[i,]$consensus_minor_allele,
                                                  "caller_agreement" = consensus_cnv_alleles_bed_df[i,]$caller_agreement,
                                                  "allele_caller_agreement" = consensus_cnv_alleles_bed_df[i,]$allele_caller_agreement))
      }
    }
  }
  return(agreement_resolved_cnv_df)
}


score_segments <- function(resolved_cnv_alleles_df) {
  
  # First calculate the size of each segment
  resolved_cnv_alleles_df <- resolved_cnv_alleles_df %>%
                              mutate(segment_size = as.numeric(end) - as.numeric(start))

  segment_confidence_score_df <- tibble()

  # For each segment multiply the segment size by the number of callers that agreed on total CN and alleles
  for(i in 1:nrow(resolved_cnv_alleles_df)) {

    consensus_cn_score_multiplier <- as.numeric(length(str_split(string = resolved_cnv_alleles_df$caller_agreement[i], 
                                                                 pattern = ",", 
                                                                 simplify = T)))
    consensus_alleles_score_multiplier <- as.numeric(length(str_split(string = resolved_cnv_alleles_df$allele_caller_agreement[i],
                                                                      pattern = ",", 
                                                                      simplify = T)))
  
    segment_confidence_score_df <- rbind(segment_confidence_score_df,
                                         tibble("consensus_cn_segment_score" = consensus_cn_score_multiplier * resolved_cnv_alleles_df$segment_size[i],
                                                "consensus_alleles_segment_score" = consensus_alleles_score_multiplier * resolved_cnv_alleles_df$segment_size[i]))
  }

  scored_resolved_cnv_alleles_df <- cbind(resolved_cnv_alleles_df, segment_confidence_score_df)
  return(scored_resolved_cnv_alleles_df)
}


chain_segments <- function(scored_resolved_cnv_alleles_df) {
  
  # Create output data frame that holds fully chained CNV segment BED lines
  chained_consensus_cnv_alleles_df <- tibble()
  
  # Loop through all entries by chromosome
  for(chromosome in unique(scored_resolved_cnv_alleles_df$chrom)) {
    
    # Subset whole CNV dataset by chromosome
    chrom_specific_subset_df <- scored_resolved_cnv_alleles_df[scored_resolved_cnv_alleles_df$chrom == chromosome,]
  
    # Create a tibble of all segments with consistent total copy number and major/minor alleles
    consistent_cn_alleles_segments_df <- tibble()
  
    for(i in 1:nrow(chrom_specific_subset_df)) {

      # Add record to start segment chaining if consistent df is empty, if not check next record
      if(is_empty(consistent_cn_alleles_segments_df)) {

        consistent_cn_alleles_segments_df <- rbind(consistent_cn_alleles_segments_df,
                                                   tibble("chrom" = chrom_specific_subset_df[i,]$chrom,
                                                      	  "start" = chrom_specific_subset_df[i,]$start,
                                                      	  "end" = chrom_specific_subset_df[i,]$end,
                                                      	  "consensus_total_cn" = chrom_specific_subset_df[i,]$consensus_total_cn,
                                                      	  "consensus_major_allele" = chrom_specific_subset_df[i,]$consensus_major_allele,
                                                      	  "consensus_minor_allele" = chrom_specific_subset_df[i,]$consensus_minor_allele,
                                                      	  "segment_size" = chrom_specific_subset_df[i,]$segment_size,
                                                      	  "consensus_cn_segment_score" = chrom_specific_subset_df[i,]$consensus_cn_segment_score,
                                                      	  "consensus_alleles_segment_score" = chrom_specific_subset_df[i,]$consensus_alleles_segment_score))
      }
    
      # Check if next record has consistent total copy number and major/minor alleles as previous record
      last_segment_index <- nrow(consistent_cn_alleles_segments_df)

      if(chrom_specific_subset_df[i,]$consensus_total_cn == consistent_cn_alleles_segments_df[last_segment_index,]$consensus_total_cn &
         chrom_specific_subset_df[i,]$consensus_major_allele == consistent_cn_alleles_segments_df[last_segment_index,]$consensus_major_allele &
         chrom_specific_subset_df[i,]$consensus_minor_allele == consistent_cn_alleles_segments_df[last_segment_index,]$consensus_minor_allele) {

        # If segment is consistent with previous, create new record that merges the segment ends
        chained_consistent_cn_alleles_segment <- tibble("chrom" = consistent_cn_alleles_segments_df[last_segment_index,]$chrom,
                                                        "start" = consistent_cn_alleles_segments_df[last_segment_index,]$start,
                                                        "end" = chrom_specific_subset_df[i,]$end,
                                                        "consensus_total_cn" = consistent_cn_alleles_segments_df[last_segment_index,]$consensus_total_cn,
                                                        "consensus_major_allele" = consistent_cn_alleles_segments_df[last_segment_index,]$consensus_major_allele,
                                                        "consensus_minor_allele" = consistent_cn_alleles_segments_df[last_segment_index,]$consensus_minor_allele,
                                                        "segment_size" = sum(consistent_cn_alleles_segments_df[last_segment_index,]$segment_size, chrom_specific_subset_df[i,]$segment_size),
                                                        "consensus_cn_segment_score" = sum(consistent_cn_alleles_segments_df[last_segment_index,]$consensus_cn_segment_score, chrom_specific_subset_df[i,]$consensus_cn_segment_score),
                                                        "consensus_alleles_segment_score" = sum(consistent_cn_alleles_segments_df[last_segment_index,]$consensus_alleles_segment_score, chrom_specific_subset_df[i,]$consensus_alleles_segment_score))

        # Replace the segment record with new chained segment
        consistent_cn_alleles_segments_df[last_segment_index,] <- chained_consistent_cn_alleles_segment
        
      } else {
      
        # If not consistent, add new segment as a record to then check for chaining
        consistent_cn_alleles_segments_df <- rbind(consistent_cn_alleles_segments_df,
                                                   tibble("chrom" = chrom_specific_subset_df[i,]$chrom,
                                                      	  "start" = chrom_specific_subset_df[i,]$start,
                                                          "end" = chrom_specific_subset_df[i,]$end,
                                                          "consensus_total_cn" = chrom_specific_subset_df[i,]$consensus_total_cn,
                                                          "consensus_major_allele" = chrom_specific_subset_df[i,]$consensus_major_allele,
                                                          "consensus_minor_allele" = chrom_specific_subset_df[i,]$consensus_minor_allele,
                                                          "segment_size" = chrom_specific_subset_df[i,]$segment_size,
                                                          "consensus_cn_segment_score" = chrom_specific_subset_df[i,]$consensus_cn_segment_score,
                                                          "consensus_alleles_segment_score" = chrom_specific_subset_df[i,]$consensus_alleles_segment_score))
      }
    }
  
    # After finding chained segments, add all records per chromosome to the final chained consensus df
    chained_consensus_cnv_alleles_df <- rbind(chained_consensus_cnv_alleles_df, consistent_cn_alleles_segments_df)
  }
  
  return(chained_consensus_cnv_alleles_df)
}


segment_classification_and_finalization <- function(chained_scored_resolved_df, sex_of_sample) {
  # First, remove any segments that are only 1bp long
  orphan_segment_indices <- c()
  final_cnv_segments_df <- tibble()

  for(i in 1:nrow(chained_scored_resolved_df)) {
  
    if(chained_scored_resolved_df[i,]$segment_size == 1) {
      orphan_segment_indices <- append(x = orphan_segment_indices,
                                       values = i)
    }
  }

  # Capture any orphaned segments and remove them from the final output segments df
  if(!is.null(orphan_segment_indices)) {

    orphaned_segments <- chained_scored_resolved_df[orphan_segment_indices,]
    final_cnv_segments_df <- setdiff(x = chained_scored_resolved_df,
                                     y = orphaned_segments)
  } else {

    final_cnv_segments_df <- chained_scored_resolved_df
  }

  # Second, generate final confidence designation based on aggregate of scores over full distance of segment
  segment_final_confidence_rating <- tibble()

  for(i in 1:nrow(final_cnv_segments_df)) {
  
    cn_segment_final_confidence <- round(final_cnv_segments_df[i,]$consensus_cn_segment_score / final_cnv_segments_df[i,]$segment_size,
                                         digits = 3)
    alleles_segment_final_confidence <- round(final_cnv_segments_df[i,]$consensus_alleles_segment_score / final_cnv_segments_df[i,]$segment_size,
                                              digits = 3)
  
    final_segment_confidence <- round((cn_segment_final_confidence + alleles_segment_final_confidence) / 2,
                                       digits = 3)
  
    # Give a final rating based on final confidence score within a range
    if(between(x = final_segment_confidence, left = 1.000, right = 1.333)) {
    
      segment_final_confidence_rating <- rbind(segment_final_confidence_rating,
                                               tibble("consensus_conf_rating" = "*"))
      
    } else if(between(x = final_segment_confidence, left = 1.334, right = 1.666)) {
    
      segment_final_confidence_rating <- rbind(segment_final_confidence_rating,
                                               tibble("consensus_conf_rating" = "**"))
      
    } else if(between(x = final_segment_confidence, left = 1.667, right = 2.000)) {
    
      segment_final_confidence_rating <- rbind(segment_final_confidence_rating,
                                               tibble("consensus_conf_rating" = "***"))
      
    } else if(final_segment_confidence >= 2.001) {
    
      segment_final_confidence_rating <- rbind(segment_final_confidence_rating,
                                               tibble("consensus_conf_rating" = "****"))
    }
  }

  final_cnv_segments_df <- cbind(final_cnv_segments_df, segment_final_confidence_rating)

  # Third, add classification to each segment based on size, total copy number, and allele status
  final_segment_classification <- tibble()

  for(i in 1:nrow(final_cnv_segments_df)) {
    segment_class <- ""
    segment_type <- ""
    segment_allele_status <- ""
  
    # Broad or focal segment?
    if(final_cnv_segments_df[i,]$segment_size < 3000000) {
      segment_class <- "focal"
    } else {
      segment_class <- "broad"
    }
  
    # Need to consider autosomes and sex chromosomes separately
    if(final_cnv_segments_df[i,]$chrom == "chrX") {
    
      # Use sample sex provided to easily determine segment type and allele status
      if(sex_of_sample == "male") {
      
        # Determine type
        if(as.numeric(final_cnv_segments_df[i,]$consensus_total_cn) == 2) {
          segment_type <- "gain"
        } else if(as.numeric(final_cnv_segments_df[i,]$consensus_total_cn) > 2) {
          segment_type <- "amplification"
        } else if(as.numeric(final_cnv_segments_df[i,]$consensus_total_cn) < 1) {
          segment_type <- "deletion"
        } else {
          segment_type <- "neutral"
        }
      
        # Determine allele status
        if(as.numeric(final_cnv_segments_df[i,]$consensus_major_allele) == 0 & as.numeric(final_cnv_segments_df[i,]$consensus_minor_allele) == 0) {
          segment_allele_status <- "hom_del"
        } else if(as.numeric(final_cnv_segments_df[i,]$consensus_major_allele) != 0 & as.numeric(final_cnv_segments_df[i,]$consensus_minor_allele) == 0) {
          segment_allele_status <- "loh"
        } else {
          segment_allele_status <- "het"
        }
      
      } else if(sex_of_sample == "female") {
      
        # Determine type
        if(as.numeric(final_cnv_segments_df[i,]$consensus_total_cn) == 3) {
          segment_type <- "gain"
        } else if(as.numeric(final_cnv_segments_df[i,]$consensus_total_cn) > 3) {
          segment_type <- "amplification"
        } else if(as.numeric(final_cnv_segments_df[i,]$consensus_total_cn) < 2) {
          segment_type <- "deletion"
        } else {
          segment_type <- "neutral"
        }
      
        # Determine allele status
        if(as.numeric(final_cnv_segments_df[i,]$consensus_major_allele) == 0 & as.numeric(final_cnv_segments_df[i,]$consensus_minor_allele) == 0) {
          segment_allele_status <- "hom_del"
        } else if(as.numeric(final_cnv_segments_df[i,]$consensus_major_allele) != 0 & as.numeric(final_cnv_segments_df[i,]$consensus_minor_allele) == 0) {
          segment_allele_status <- "loh"
        } else {
          segment_allele_status <- "het"
        }
      }
    
    } else if(final_cnv_segments_df[i,]$chrom == "chrY") {
    
      # Use sample sex provided to easily determine segment type and allele status
      if(sex_of_sample == "male") {
      
        # Determine type
        if(as.numeric(final_cnv_segments_df[i,]$consensus_total_cn) == 2) {
          segment_type <- "gain"
        } else if(as.numeric(final_cnv_segments_df[i,]$consensus_total_cn) > 2) {
          segment_type <- "amplification"
        } else if(as.numeric(final_cnv_segments_df[i,]$consensus_total_cn) < 1) {
          segment_type <- "deletion"
        } else {
          segment_type <- "neutral"
        }
      
        # Determine allele status
        if(as.numeric(final_cnv_segments_df[i,]$consensus_major_allele) == 0 & as.numeric(final_cnv_segments_df[i,]$consensus_minor_allele) == 0) {
          segment_allele_status <- "hom_del"
        } else if(as.numeric(final_cnv_segments_df[i,]$consensus_major_allele) != 0 & as.numeric(final_cnv_segments_df[i,]$consensus_minor_allele) == 0) {
          segment_allele_status <- "loh"
        } else {
          segment_allele_status <- "het"
        }
      
      } else if(sex_of_sample == "female") {
      
        # Determine type
        if(as.numeric(final_cnv_segments_df[i,]$consensus_total_cn) == 1) {
          segment_type <- "gain"
        } else if(as.numeric(final_cnv_segments_df[i,]$consensus_total_cn) > 1) {
          segment_type <- "amplification"
        } else {
          segment_type <- "neutral"
        }
      
        # Determine allele status
        if(as.numeric(final_cnv_segments_df[i,]$consensus_major_allele) == 0 & as.numeric(final_cnv_segments_df[i,]$consensus_minor_allele) == 0) {
          segment_allele_status <- "hom_del"
        } else if(as.numeric(final_cnv_segments_df[i,]$consensus_major_allele) != 0 & as.numeric(final_cnv_segments_df[i,]$consensus_minor_allele) == 0) {
          segment_allele_status <- "loh"
        } else {
          segment_allele_status <- "het"
        }
      }
    
    } else {
    
      # Determine type
      if(as.numeric(final_cnv_segments_df[i,]$consensus_total_cn) == 3) {
        segment_type <- "gain"
      } else if(as.numeric(final_cnv_segments_df[i,]$consensus_total_cn) > 3) {
        segment_type <- "amplification"
      } else if(as.numeric(final_cnv_segments_df[i,]$consensus_total_cn) < 2) {
        segment_type <- "deletion"
      } else {
        segment_type <- "neutral"
      }
    
      # Determine allele status
      if(as.numeric(final_cnv_segments_df[i,]$consensus_major_allele) == 0 & as.numeric(final_cnv_segments_df[i,]$consensus_minor_allele) == 0) {
        segment_allele_status <- "hom_del"
      } else if(as.numeric(final_cnv_segments_df[i,]$consensus_major_allele) != 0 & as.numeric(final_cnv_segments_df[i,]$consensus_minor_allele) == 0) {
        segment_allele_status <- "loh"
      } else {
        segment_allele_status <- "het"
      }
    }
  
    final_segment_classification <- rbind(final_segment_classification,
                                        tibble("type" = segment_type,
                                               "allele_status" = segment_allele_status,
                                               "class" = segment_class))
  }

  final_cnv_segments_df <- cbind(final_cnv_segments_df, final_segment_classification)

  # Finally, isolate desired columns and output the final df
  final_high_quality_cnv_bed <- final_cnv_segments_df %>%
                                select(chrom, start, end,
                                       consensus_total_cn, consensus_major_allele, consensus_minor_allele,
                                       type, class, allele_status, consensus_conf_rating)
  
  return(final_high_quality_cnv_bed)
}

#########################
#####   Execution   #####

# Accept command line arguments as input
input_args <- commandArgs(trailingOnly = T)

sample_sex <- input_args[2]
output_bed_filename <- input_args[3]
mechansim <- input_args[4]

print("Reading input consensus merged CNV and allele BED file.....")
consensus_cnv_alleles_bed_file <- read_delim(file = input_args[1],
                                            delim = "\t",
                                            col_names = T,
                                            col_types = "cccccccccccccc")

print("Resolving all segments with no agreement between callers.....")
resolved_cnv_alleles <- segment_disagreement_resolution(consensus_cnv_alleles_bed_df = consensus_cnv_alleles_bed_file, consensus_mechanism = mechansim)

print("Scoring each segment based on length and caller agreement.....")
scored_resolved_cnv_alleles <- score_segments(resolved_cnv_alleles_df = resolved_cnv_alleles)

print("Chaining all segments that have consistent total CN and major/minor alleles.....")
chained_scored_resolved_cnv_alleles <- chain_segments(scored_resolved_cnv_alleles_df = scored_resolved_cnv_alleles)

print("Classifying and finalizing all segments based on confidence, size, type, and allele status.....")
final_high_quality_cnv_output <- segment_classification_and_finalization(chained_scored_resolved_df = chained_scored_resolved_cnv_alleles, sex_of_sample = sample_sex)

write.table(x = final_high_quality_cnv_output,
            file = output_bed_filename,
            sep = "\t",
            col.names = TRUE,
            quote = FALSE,
            row.names = FALSE)
