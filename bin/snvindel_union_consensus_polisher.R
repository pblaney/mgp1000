#!/usr/bin/env Rscript

# This script accepts set of SNV or InDel VCF files as input and generates
# a union consensus mutation file

#########################
#####   Libraries   #####

suppressPackageStartupMessages(library(devgru))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(VariantAnnotation))
suppressPackageStartupMessages(library(gUtils))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doParallel))

Sys.setenv(DEFAULT_BSGENOME = 'BSgenome.Hsapiens.UCSC.hg38::Hsapiens')
options(scipen = 999)

#########################
#####   Functions   #####

# Explicitly calculate the VAF for mutations of various flavors (strelka, svaba)
calc_vaf <- function(mut_record, caller, mut_type) {
  
  if(caller == "strelka" & mut_type == "snv") {
    # Strelka SNVs
    # Does not directly report VAF per mutation, instead it reports the depth DP ["Read depth for tier1 (used+filtered)"]
    # also reports per nucleotide tier read support AU, CU, TU, GU ["Number of 'A/C/G/T' alleles used in tiers 1,2"]
    # The VAF will be calculated as: [AU_TIER1, CU_TIER1, TU_TIER1, GU_TIER1] / DP
    
    # TODO: need to also account for non-union/default format Strelka DT/GR
    alt_read_depth <- as.data.frame(mut_record)[,colnames(mut_record) == paste0("STRELKA_", mut_record$ALT, "U_TIER1_TUMOR")]

    mut_vaf <- alt_read_depth / mut_record$STRELKA_DP_TUMOR
    
  } else if(caller == "strelka" & mut_type == "indel") {
    # Strelka InDels
    # Like SNVs, reports depth DP and indel tier read support TIR ["Reads strongly supporting indel allele for tiers 1,2"]
    # The VAF will be calculated as: TIR_TIER1 / DP
    mut_vaf <- mut_record$STRELKA_TIR_TIER1_TUMOR / mut_record$STRELKA_DP_TUMOR
    
  } else if(caller == "svaba" & mut_type == "indel") {
    # SvABA InDels
    # Does not directly report VAF per mutation, instead it reports depth DP ["Depth of coverage: Number of reads covering site."]
    # also reports allele depth AD ["Allele depth: Number of reads supporting the variant"]
    # The VAF will be calculated as: AD / DP
    mut_vaf <- mut_record$SVABA_AD_TUMOR / mut_record$SVABA_DP_TUMOR
  }
  
  # Output the VAF for the mutation record
  return(mut_vaf)
}


#########################
#####   Execution   #####

# Accept command line arguments as input
input_args <- commandArgs(trailingOnly = T)

vcf_file_dir <- input_args[1]

mut_type <- input_args[2]

output_dir <- input_args[3]

gene_gtf_file <- input_args[4]

threads <- input_args[5]

# Set number of threads for foreach loop and file writing
registerDoParallel(cores = threads)

# Find all input files in directory provided
message("\nScanning path ", vcf_file_dir, " for input ", stringr::str_to_upper(mut_type), " VCFs.....")
mutect_vcf_input <- list.files(path = vcf_file_dir,
                               pattern = paste0("*.mutect.somatic.", mut_type, ".vcf.gz"))

strelka_vcf_input <- list.files(path = vcf_file_dir,
                                pattern = paste0("*.strelka.somatic.", mut_type, ".vcf.gz"))

varscan_vcf_input <- list.files(path = vcf_file_dir,
                                pattern = paste0("*.varscan.somatic.", mut_type, ".vcf.gz"))

# Edge case: include SvABA for InDel consensus
if(mut_type == "indel") {
  svaba_vcf_input <- list.files(path = vcf_file_dir,
                                pattern = "*.svaba.somatic.indel.vcf.gz")
}

# Get sample names and gather input file pairs
message("Collecting possible per sample VCFs.....")
mutect_samples <- stringr::str_remove(string = mutect_vcf_input, pattern = ".mutect.somatic.*.vcf.gz") %>%
  tibble::as_tibble_col(column_name = "sample")
message(nrow(mutect_samples), " Mutect VCF(s) detected....")

strelka_samples <- stringr::str_remove(string = strelka_vcf_input, pattern = ".strelka.somatic.*.vcf.gz") %>%
  tibble::as_tibble_col(column_name = "sample")
message(nrow(strelka_samples), " Strelka VCF(s) detected....")

varscan_samples <- stringr::str_remove(string = varscan_vcf_input, pattern = ".varscan.somatic.*.vcf.gz") %>%
  tibble::as_tibble_col(column_name = "sample")
message(nrow(varscan_samples), " Varscan VCF(s) detected....")

if(mut_type == "snv") {
  sample_set <- dplyr::inner_join(x = mutect_samples, y = strelka_samples, by = "sample") %>%
    dplyr::inner_join(y = varscan_samples, by = "sample")
  message("\n", nrow(sample_set), " samples with triplet VCF(s)....")
  
# Edge case: include SvABA for InDel consensus
} else if(mut_type == "indel") {
  svaba_samples <- stringr::str_remove(string = svaba_vcf_input, pattern = ".svaba.somatic.indel.vcf.gz") %>%
    tibble::as_tibble_col(column_name = "sample")
  message(nrow(svaba_samples), " SvABA VCF(s) detected....")
  
  sample_set <- dplyr::inner_join(x = mutect_samples, y = strelka_samples, by = "sample") %>%
    dplyr::inner_join(y = varscan_samples, by = "sample") %>%
    dplyr::inner_join(y = svaba_samples, by = "sample")
  message("\n", nrow(sample_set), " samples with quadra VCF(s)....")
}

# Read in genes
genes <- get_genes_shortcut(gtf_file_path = gene_gtf_file)

# Loop through all per sample triplet VCFs to create consensus mutation set
message("Running union consensus using ", threads, " threads.....")
for(i in 1:nrow(sample_set)) {
  
  # Get paths to VCFs
  mutect_vcf <- mutect_vcf_input[stringr::str_detect(string = mutect_vcf_input, pattern = sample_set$sample[i]) %>% which()]
  strelka_vcf <- strelka_vcf_input[stringr::str_detect(string = strelka_vcf_input, pattern = sample_set$sample[i]) %>% which()]
  varscan_vcf <- varscan_vcf_input[stringr::str_detect(string = varscan_vcf_input, pattern = sample_set$sample[i]) %>% which()]
  
  # Edge case: include SvABA for InDel consensus
  if(mut_type == "indel") {
    svaba_vcf <- svaba_vcf_input[stringr::str_detect(string = svaba_vcf_input, pattern = sample_set$sample[i]) %>% which()]
  }
  
  # Read in a VCF file and convert to GRanges object
  mutect_gr <- read_vcf_file(vcf_file_path = paste0(vcf_file_dir, mutect_vcf),
                             tumor_sample = stringr::str_split(string = sample_set$sample[i], pattern = "_vs_", simplify = T)[,1],
                             normal_sample = stringr::str_split(string = sample_set$sample[i], pattern = "_vs_", simplify = T)[,2],
                             caller = "mutect",
                             mut_type = mut_type)
  # Add prefix to column names to keep uniqueness
  colnames(mcols(mutect_gr))[8:length(mcols(mutect_gr))] <- stringr::str_c("MUTECT_",
                                                                           colnames(mcols(mutect_gr))[8:length(mcols(mutect_gr))])
  
  strelka_gr <- read_vcf_file(vcf_file_path = paste0(vcf_file_dir, strelka_vcf),
                              tumor_sample = stringr::str_split(string = sample_set$sample[i], pattern = "_vs_", simplify = T)[,1],
                              normal_sample = stringr::str_split(string = sample_set$sample[i], pattern = "_vs_", simplify = T)[,2],
                              caller = "strelka",
                              mut_type = mut_type)
  # Add prefix to column names to keep uniqueness
  colnames(mcols(strelka_gr))[8:length(mcols(strelka_gr))] <- stringr::str_c("STRELKA_",
                                                                             colnames(mcols(strelka_gr))[8:length(mcols(strelka_gr))])
  
  varscan_gr <- read_vcf_file(vcf_file_path = paste0(vcf_file_dir, varscan_vcf),
                              tumor_sample = stringr::str_split(string = sample_set$sample[i], pattern = "_vs_", simplify = T)[,1],
                              normal_sample = stringr::str_split(string = sample_set$sample[i], pattern = "_vs_", simplify = T)[,2],
                              caller = "varscan",
                              mut_type = mut_type)
  # Add prefix to column names to keep uniqueness
  colnames(mcols(varscan_gr))[8:length(mcols(varscan_gr))] <- stringr::str_c("VARSCAN_",
                                                                             colnames(mcols(varscan_gr))[8:length(mcols(varscan_gr))])
  
  # Edge case: include SvABA for InDel consensus
  if(mut_type == "indel") {
    svaba_gr <- read_vcf_file(vcf_file_path = paste0(vcf_file_dir, svaba_vcf),
                              tumor_sample = stringr::str_split(string = sample_set$sample[i], pattern = "_vs_", simplify = T)[,1],
                              normal_sample = stringr::str_split(string = sample_set$sample[i], pattern = "_vs_", simplify = T)[,2],
                              caller = "svaba",
                              mut_type = "indel")
    # Add prefix to column names to keep uniqueness
    colnames(mcols(svaba_gr))[8:length(mcols(svaba_gr))] <- stringr::str_c("SVABA_",
                                                                           colnames(mcols(svaba_gr))[8:length(mcols(svaba_gr))])
  }
  
  # Merge all the calls into single, sorted GR obj
  if(mut_type == "snv") {
    union_gr <- gUtils::grbind(mutect_gr, strelka_gr, varscan_gr)
  } else if(mut_type == "indel") {
    union_gr <- gUtils::grbind(mutect_gr, strelka_gr, varscan_gr, svaba_gr)
  }
  union_gr <- GenomicRanges::sort.GenomicRanges(union_gr)
  
  # Refactor seqinfo for consistency
  union_gr <- gr_refactor_seqs(input_gr = union_gr)
  
  # Loop through all merged calls to find and merge consensus records (i.e. same call from multiple callers)
  # Split each input loop by chromosome
  chrom_iter_list <- as.character(union_gr@seqnames@values)
  
  final_union_consensus_gr <- foreach(x = 1:length(chrom_iter_list), .combine = grbind, .packages = "gUtils") %dopar% {
    
    union_per_chrom_gr <- union_gr %Q% (seqnames == chrom_iter_list[x])
    
    # Create empty GRanges
    union_consensus_gr <- GenomicRanges::GRanges()
    
    for(j in 1:length(union_per_chrom_gr)) {
      
      # Check each record for any overlapping records by matching chrom, start, end, patient, ref, and alt columns
      test_for_overlap <- gUtils::gr.findoverlaps(query = union_per_chrom_gr[j],
                                                  subject = union_per_chrom_gr,
                                                  scol = colnames(mcols(union_per_chrom_gr)),
                                                  by = c("PATIENT", "REF", "ALT"))
      
      # overlap test GR will be 1,2, or 3
      if(length(test_for_overlap) > 1) {
        # set new base output that contains all columns and those specific to first caller in consensus
        consensus_gr <- test_for_overlap[1]
        
        # concat caller strings for new caller column
        consensus_gr$CALLER <- stringr::str_flatten(sort(test_for_overlap$CALLER), collapse = ",")
        
        # now cycle through the other records to grab caller specific metrics
        for(k in 2:length(test_for_overlap)) {
          # which caller is the record
          which_caller <- stringr::str_to_upper(test_for_overlap$CALLER[k])
          
          # which GR VCF metadata column indices correspond to that caller
          caller_col_idx <- grep(pattern = which_caller, x = colnames(mcols(test_for_overlap)))
          
          # replace the placeholder NA's with caller specific metadata
          mcols(consensus_gr)[caller_col_idx] <- mcols(test_for_overlap[k, caller_col_idx])
        }
        
        # Add the consensus record to the output, skip if already in there
        if(sum(gUtils::gr.in(query = union_consensus_gr, subject = consensus_gr)) > 0) {
          next
        } else {
          union_consensus_gr <- gUtils::grbind(union_consensus_gr, consensus_gr)
        }
        
        # Single caller situations
      } else if(length(test_for_overlap) == 1) {
        # Add the singleton record to the output
        singleton_gr <- test_for_overlap[1]
        union_consensus_gr <- gUtils::grbind(union_consensus_gr, singleton_gr)
      }
    }
    
    # Return output of the foreach loops, each GR obj will be concatenated
    union_consensus_gr
  }
  
  # Annotated each mutation with nearest gene for rapid identification of potential driver muts
  message("Annotating mutations by nearest gene.....")
  final_union_consensus_gr$nearest_gene <- genes$gene_name[IRanges::nearest(x = final_union_consensus_gr,
                                                                            subject = gUtils::gr.stripstrand(genes))]
  
  # Slim down object before output (rm query.id, subject.id)
  final_union_consensus_gr <- final_union_consensus_gr[,c(-1,-2)]

  # Convert GRanges to DT
  final_union_consensus_dt <- gUtils::gr2dt(x = final_union_consensus_gr)

  # Now loop through all union consensus calls to report VAFs in single column
  final_union_consensus_vaf_metrics_dt <- foreach(i = 1:length(chrom_iter_list), .combine = rrbind, .packages = "gUtils") %dopar% {

    # Split each input loop by chromosome
    final_union_consensus_per_chrom_dt <- final_union_consensus_dt[final_union_consensus_dt$seqnames == chrom_iter_list[i]]

    vaf_info <- data.table(vaf_combo = rep(NA, nrow(final_union_consensus_per_chrom_dt)),
                           vaf_mean = rep(NA, nrow(final_union_consensus_per_chrom_dt)))

    for(j in 1:nrow(final_union_consensus_per_chrom_dt)) {
      
      # First, determine which callers found the variant to start the decision tree
      caller_string <- final_union_consensus_per_chrom_dt[j, CALLER]

      # Use the string of callers to find correct VAF combo and mean
      if(caller_string == "mutect,strelka,svaba,varscan") {
        # 4 CALLER CONSENSUS
        mutect_vaf <- round(final_union_consensus_per_chrom_dt[j, MUTECT_VAF],
                            digits = 4)
        strelka_vaf <- round(calc_vaf(mut_record = final_union_consensus_per_chrom_dt[j,], caller = "strelka", mut_type = mut_type),
                             digits = 4)
        varscan_vaf <- round(as.numeric(str_remove(string = final_union_consensus_per_chrom_dt[j, VARSCAN_FREQ_TUMOR], pattern = "%")) / 100,
                             digits = 4)
        svaba_vaf <- round(calc_vaf(mut_record = final_union_consensus_per_chrom_dt[j,], caller = "svaba", mut_type = mut_type),
                           digits = 4)
        
        # Ouput string of individual VAFs in same order as CALLER and a mean of the values
        vaf_combo_string <- str_c(mutect_vaf, strelka_vaf, svaba_vaf, varscan_vaf,
                                  sep = ",")
        mean_vaf <- round(mean(c(mutect_vaf, strelka_vaf, svaba_vaf, varscan_vaf)),
                          digits = 4)
        
      } else if(caller_string == "mutect,strelka,svaba") {
        # 3 CALLER CONSENSUS
        mutect_vaf <- round(final_union_consensus_per_chrom_dt[j, MUTECT_VAF],
                            digits = 4)
        strelka_vaf <- round(calc_vaf(mut_record = final_union_consensus_per_chrom_dt[j,], caller = "strelka", mut_type = mut_type),
                             digits = 4)
        svaba_vaf <- round(calc_vaf(mut_record = final_union_consensus_per_chrom_dt[j,], caller = "svaba", mut_type = mut_type),
                           digits = 4)
        
        # Ouput string of individual VAFs in same order as CALLER and a mean of the values
        vaf_combo_string <- str_c(mutect_vaf, strelka_vaf, svaba_vaf, 
                                  sep = ",")
        mean_vaf <- round(mean(c(mutect_vaf, strelka_vaf, svaba_vaf)),
                          digits = 4)
        
      } else if(caller_string == "mutect,strelka,varscan") {
        # 3 CALLER CONSENSUS
        mutect_vaf <- round(final_union_consensus_per_chrom_dt[j, MUTECT_VAF],
                            digits = 4)
        strelka_vaf <- round(calc_vaf(mut_record = final_union_consensus_per_chrom_dt[j,], caller = "strelka", mut_type = mut_type),
                             digits = 4)
        varscan_vaf <- round(as.numeric(stringr::str_remove(string = final_union_consensus_per_chrom_dt[j, VARSCAN_FREQ_TUMOR], pattern = "%")) / 100,
                             digits = 4)
        
        # Ouput string of individual VAFs in same order as CALLER and a mean of the values
        vaf_combo_string <- stringr::str_c(mutect_vaf, strelka_vaf, varscan_vaf, 
                                           sep = ",")
        mean_vaf <- round(mean(c(mutect_vaf, strelka_vaf, varscan_vaf)),
                          digits = 4)
        
      } else if(caller_string == "mutect,svaba,varscan") {
        # 3 CALLER CONSENSUS
        mutect_vaf <- round(final_union_consensus_per_chrom_dt[j, MUTECT_VAF],
                            digits = 4)
        svaba_vaf <- round(calc_vaf(mut_record = final_union_consensus_per_chrom_dt[j,], caller = "svaba", mut_type = mut_type),
                           digits = 4)
        varscan_vaf <- round(as.numeric(stringr::str_remove(string = final_union_consensus_per_chrom_dt[j, VARSCAN_FREQ_TUMOR], pattern = "%")) / 100,
                             digits = 4)
        
        # Ouput string of individual VAFs in same order as CALLER and a mean of the values
        vaf_combo_string <- stringr::str_c(mutect_vaf, svaba_vaf, varscan_vaf,
                                           sep = ",")
        mean_vaf <- round(mean(c(mutect_vaf, svaba_vaf, varscan_vaf)),
                          digits = 4)
        
      } else if(caller_string == "strelka,svaba,varscan") {
        # 3 CALLER CONSENSUS
        strelka_vaf <- round(calc_vaf(mut_record = final_union_consensus_per_chrom_dt[j,], caller = "strelka", mut_type = mut_type),
                             digits = 4)
        svaba_vaf <- round(calc_vaf(mut_record = final_union_consensus_per_chrom_dt[j,], caller = "svaba", mut_type = mut_type),
                           digits = 4)
        varscan_vaf <- round(as.numeric(stringr::str_remove(string = final_union_consensus_per_chrom_dt[j, VARSCAN_FREQ_TUMOR], pattern = "%")) / 100,
                             digits = 4)
        
        # Ouput string of individual VAFs in same order as CALLER and a mean of the values
        vaf_combo_string <- stringr::str_c(strelka_vaf, svaba_vaf, varscan_vaf,
                                           sep = ",")
        mean_vaf <- round(mean(c(strelka_vaf, svaba_vaf, varscan_vaf)),
                          digits = 4)
        
      } else if(caller_string == "mutect,strelka") {
        # 2 CALLER CONSENSUS
        mutect_vaf <- round(final_union_consensus_per_chrom_dt[j, MUTECT_VAF],
                            digits = 4)
        strelka_vaf <- round(calc_vaf(mut_record = final_union_consensus_per_chrom_dt[j,], caller = "strelka", mut_type = mut_type),
                             digits = 4)
        
        # Ouput string of individual VAFs in same order as CALLER and a mean of the values
        vaf_combo_string <- stringr::str_c(mutect_vaf, strelka_vaf,
                                           sep = ",")
        mean_vaf <- round(mean(c(mutect_vaf, strelka_vaf)),
                          digits = 4)
        
      } else if(caller_string == "mutect,svaba") {
        # 2 CALLER CONSENSUS
        mutect_vaf <- round(final_union_consensus_per_chrom_dt[j, MUTECT_VAF],
                            digits = 4)
        svaba_vaf <- round(calc_vaf(mut_record = final_union_consensus_per_chrom_dt[j,], caller = "svaba", mut_type = mut_type),
                           digits = 4)
        
        # Ouput string of individual VAFs in same order as CALLER and a mean of the values
        vaf_combo_string <- stringr::str_c(mutect_vaf, svaba_vaf,
                                           sep = ",")
        mean_vaf <- round(mean(c(mutect_vaf, svaba_vaf)),
                          digits = 4)
        
      } else if(caller_string == "mutect,varscan") {
        # 2 CALLER CONSENSUS
        mutect_vaf <- round(final_union_consensus_per_chrom_dt[j, MUTECT_VAF],
                            digits = 4)
        varscan_vaf <- round(as.numeric(stringr::str_remove(string = final_union_consensus_per_chrom_dt[j, VARSCAN_FREQ_TUMOR], pattern = "%")) / 100,
                             digits = 4)
        
        # Ouput string of individual VAFs in same order as CALLER and a mean of the values
        vaf_combo_string <- stringr::str_c(mutect_vaf, varscan_vaf,
                                           sep = ",")
        mean_vaf <- round(mean(c(mutect_vaf, varscan_vaf)),
                          digits = 4)
        
      } else if(caller_string == "strelka,svaba") {
        # 2 CALLER CONSENSUS
        strelka_vaf <- round(calc_vaf(mut_record = final_union_consensus_per_chrom_dt[j,], caller = "strelka", mut_type = mut_type),
                             digits = 4)
        svaba_vaf <- round(calc_vaf(mut_record = final_union_consensus_per_chrom_dt[j,], caller = "svaba", mut_type = mut_type),
                           digits = 4)
        
        # Ouput string of individual VAFs in same order as CALLER and a mean of the values
        vaf_combo_string <- stringr::str_c(strelka_vaf, svaba_vaf,
                                           sep = ",")
        mean_vaf <- round(mean(c(strelka_vaf, svaba_vaf)),
                          digits = 4)
        
      } else if(caller_string == "strelka,varscan") {
        # 2 CALLER CONSENSUS
        strelka_vaf <- round(calc_vaf(mut_record = final_union_consensus_per_chrom_dt[j,], caller = "strelka", mut_type = mut_type),
                             digits = 4)
        varscan_vaf <- round(as.numeric(stringr::str_remove(string = final_union_consensus_per_chrom_dt[j, VARSCAN_FREQ_TUMOR], pattern = "%")) / 100,
                             digits = 4)
        
        # Ouput string of individual VAFs in same order as CALLER and a mean of the values
        vaf_combo_string <- stringr::str_c(strelka_vaf, varscan_vaf, 
                                           sep = ",")
        mean_vaf <- round(mean(c(strelka_vaf, varscan_vaf)),
                          digits = 4)
        
      } else if(caller_string == "svaba,varscan") {
        # 2 CALLER CONSENSUS
        svaba_vaf <- round(calc_vaf(mut_record = final_union_consensus_per_chrom_dt[j,], caller = "svaba", mut_type = mut_type),
                           digits = 4)
        varscan_vaf <- round(as.numeric(stringr::str_remove(string = final_union_consensus_per_chrom_dt[j, VARSCAN_FREQ_TUMOR], pattern = "%")) / 100,
                             digits = 4)
        
        # Ouput string of individual VAFs in same order as CALLER and a mean of the values
        vaf_combo_string <- stringr::str_c(svaba_vaf, varscan_vaf,
                                           sep = ",")
        mean_vaf <- round(mean(c(svaba_vaf, varscan_vaf)),
                          digits = 4)

      } else if(caller_string == "mutect") {
        # 1 CALLER SINGLETON
        mutect_vaf <- round(final_union_consensus_per_chrom_dt[j, MUTECT_VAF],
                            digits = 4)

        # Ouput string of individual VAFs in same order as CALLER and a mean of the values
        vaf_combo_string <- mutect_vaf
        mean_vaf <- mutect_vaf

      } else if(caller_string == "strelka") {
        # 1 CALLER SINGLETON
        strelka_vaf <- round(calc_vaf(mut_record = final_union_consensus_per_chrom_dt[j,], caller = "strelka", mut_type = mut_type),
                             digits = 4)

        # Ouput string of individual VAFs in same order as CALLER and a mean of the values
        vaf_combo_string <- strelka_vaf
        mean_vaf <- strelka_vaf

      } else if(caller_string == "svaba") {
        # 1 CALLER SINGLETON
        svaba_vaf <- round(calc_vaf(mut_record = final_union_consensus_per_chrom_dt[j,], caller = "svaba", mut_type = mut_type),
                           digits = 4)

        # Ouput string of individual VAFs in same order as CALLER and a mean of the values
        vaf_combo_string <- svaba_vaf
        mean_vaf <- svaba_vaf

      } else if(caller_string == "varscan") {
        # 1 CALLER SINGLETON
        varscan_vaf <- round(as.numeric(stringr::str_remove(string = final_union_consensus_per_chrom_dt[j, VARSCAN_FREQ_TUMOR], pattern = "%")) / 100,
                             digits = 4)

        # Ouput string of individual VAFs in same order as CALLER and a mean of the values
        vaf_combo_string <- varscan_vaf
        mean_vaf <- varscan_vaf
      }
      
      # Add all VAF info into output DT record
      vaf_info$vaf_combo[j] <- vaf_combo_string
      vaf_info$vaf_mean[j] <- mean_vaf
    }
    
    # print foreach output to conclude parallel run
    vaf_info
  }

  # bind the VAF metrics with the main mutation DT
  final_union_consensus_vaf_dt <- cbind(final_union_consensus_dt, final_union_consensus_vaf_metrics_dt)

  # Write the output
  message("Writing output .txt files to ", output_dir, ".....")
  data.table::fwrite(x = final_union_consensus_dt,
                     file = paste0(output_dir, "/", sample_set$sample[i], ".hq.union.consensus.somatic.", mut_type, ".txt.gz"),
                     quote = FALSE,
                     sep = "\t",
                     row.names = FALSE,
                     col.names = TRUE,
                     logical01 = FALSE,
                     na = "NA",
                     nThread = threads)
  
  message(sample_set$sample[i], " ", stringr::str_to_upper(mut_type), "s..... D O N E\n")
}
