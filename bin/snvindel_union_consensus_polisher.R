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
message("\nSetting parallel cores for execution to ", threads, " ...")
registerDoParallel(cores = threads)

# Find all input files in directory provided
message("\nScanning path ", vcf_file_dir, " for input ", stringr::str_to_upper(mut_type), " VCFs ...")
mutect_vcf_input <- list.files(path = vcf_file_dir,
                               pattern = paste0("*.mutect.somatic.", mut_type, ".vcf.gz"))

strelka_vcf_input <- list.files(path = vcf_file_dir,
                                pattern = paste0("*.strelka.somatic.", mut_type, ".vcf.gz"))

varscan_vcf_input <- list.files(path = vcf_file_dir,
                                pattern = paste0("*.varscan.somatic.", mut_type, ".vcf.gz"))

caveman_vcf_input <- list.files(path = vcf_file_dir,
                                pattern = "*.caveman.somatic.snv.vcf.gz")

svaba_vcf_input <- list.files(path = vcf_file_dir,
                              pattern = "*.svaba.somatic.indel.vcf.gz")

# Get sample names and gather input file pairs
message("Collecting possible per sample VCFs.....")
mutect_samples <- stringr::str_remove(string = mutect_vcf_input, pattern = ".mutect.somatic.*.vcf.gz") %>%
  tibble::as_tibble_col(column_name = "sample")
message(nrow(mutect_samples), " Mutect VCF(s) detected ...")

strelka_samples <- stringr::str_remove(string = strelka_vcf_input, pattern = ".strelka.somatic.*.vcf.gz") %>%
  tibble::as_tibble_col(column_name = "sample")
message(nrow(strelka_samples), " Strelka VCF(s) detected ...")

varscan_samples <- stringr::str_remove(string = varscan_vcf_input, pattern = ".varscan.somatic.*.vcf.gz") %>%
  tibble::as_tibble_col(column_name = "sample")
message(nrow(varscan_samples), " Varscan VCF(s) detected ...")

if(mut_type == "snv" & length(caveman_vcf_input) == 0) {
  sample_set <- dplyr::inner_join(x = mutect_samples, y = strelka_samples, by = "sample") %>%
    dplyr::inner_join(y = varscan_samples, by = "sample")
  message("\n", nrow(sample_set), " samples with triplet VCF(s) ...")

  # Edge case: include CaVEMan for SNV consensus
} else if(mut_type == "snv" & length(caveman_vcf_input) == 1) {
  caveman_samples <- stringr::str_remove(string = caveman_vcf_input, pattern = ".caveman.somatic.snv.vcf.gz") %>%
    tibble::as_tibble_col(column_name = "sample")
  message(nrow(caveman_samples), " CaVEMan VCF(s) detected ...")
  
  sample_set <- dplyr::inner_join(x = mutect_samples, y = strelka_samples, by = "sample") %>%
    dplyr::inner_join(y = varscan_samples, by = "sample") %>%
    dplyr::inner_join(y = caveman_samples, by = "sample")
  message("\n", nrow(sample_set), " samples with quadra VCF(s) ...")
  
  # Edge case: include SvABA for InDel consensus
} else if(mut_type == "indel") {
  svaba_samples <- stringr::str_remove(string = svaba_vcf_input, pattern = ".svaba.somatic.indel.vcf.gz") %>%
    tibble::as_tibble_col(column_name = "sample")
  message(nrow(svaba_samples), " SvABA VCF(s) detected ...")
  
  sample_set <- dplyr::inner_join(x = mutect_samples, y = strelka_samples, by = "sample") %>%
    dplyr::inner_join(y = varscan_samples, by = "sample") %>%
    dplyr::inner_join(y = svaba_samples, by = "sample")
  message("\n", nrow(sample_set), " samples with quadra VCF(s) ...")
}

# Read in genes
genes <- get_genes_shortcut(gtf_file_path = gene_gtf_file)

# Loop through all per sample triplet VCFs to create consensus mutation set
for(i in 1:nrow(sample_set)) {
  
  # Get paths to VCFs
  mutect_vcf <- mutect_vcf_input[stringr::str_detect(string = mutect_vcf_input, pattern = sample_set$sample[i]) %>% which()]
  strelka_vcf <- strelka_vcf_input[stringr::str_detect(string = strelka_vcf_input, pattern = sample_set$sample[i]) %>% which()]
  varscan_vcf <- varscan_vcf_input[stringr::str_detect(string = varscan_vcf_input, pattern = sample_set$sample[i]) %>% which()]
  
  # Edge case: include CaVEMan for SNV consensus
  if(mut_type == "snv" & length(caveman_vcf_input) == 1) {
    caveman_vcf <- caveman_vcf_input[stringr::str_detect(string = caveman_vcf_input, pattern = sample_set$sample[i]) %>% which()]
  
    # Edge case: include SvABA for InDel consensus
  } else if(mut_type == "indel") {
    svaba_vcf <- svaba_vcf_input[stringr::str_detect(string = svaba_vcf_input, pattern = sample_set$sample[i]) %>% which()]
  }
  
  # Read in a VCF file and convert to GRanges object
  mutect_gr <- read_vcf_file(vcf_file_path = paste0(vcf_file_dir, mutect_vcf),
                             tumor_sample = stringr::str_split(string = sample_set$sample[i], pattern = "_vs_", simplify = T)[,1],
                             normal_sample = stringr::str_split(string = sample_set$sample[i], pattern = "_vs_", simplify = T)[,2],
                             caller = "mutect",
                             mut_type = mut_type)
  # Get the caller specific alternate read depth, total read depth, and VAF into common format
  mutect_vaf_dt <- get_vaf(vcf_obj = mutect_gr,
                           caller = "mutect",
                           mut_type = mut_type)
  # Add the common metadata to the GR obj
  S4Vectors::mcols(mutect_gr) <- c(S4Vectors::mcols(mutect_gr), mutect_vaf_dt)
  
  # Add prefix to column names to keep uniqueness
  colnames(mcols(mutect_gr))[8:length(mcols(mutect_gr))] <- stringr::str_c("MUTECT_",
                                                                           colnames(mcols(mutect_gr))[8:length(mcols(mutect_gr))])
  
  strelka_gr <- read_vcf_file(vcf_file_path = paste0(vcf_file_dir, strelka_vcf),
                              tumor_sample = stringr::str_split(string = sample_set$sample[i], pattern = "_vs_", simplify = T)[,1],
                              normal_sample = stringr::str_split(string = sample_set$sample[i], pattern = "_vs_", simplify = T)[,2],
                              caller = "strelka",
                              mut_type = mut_type)
  # Get the caller specific alternate read depth, total read depth, and VAF into common format
  strelka_vaf_dt <- get_vaf(vcf_obj = strelka_gr,
                           caller = "strelka",
                           mut_type = mut_type)
   # Add the common metadata to the GR obj
  S4Vectors::mcols(strelka_gr) <- c(S4Vectors::mcols(strelka_gr), strelka_vaf_dt)
  # Add prefix to column names to keep uniqueness
  colnames(mcols(strelka_gr))[8:length(mcols(strelka_gr))] <- stringr::str_c("STRELKA_",
                                                                             colnames(mcols(strelka_gr))[8:length(mcols(strelka_gr))])
  
  varscan_gr <- read_vcf_file(vcf_file_path = paste0(vcf_file_dir, varscan_vcf),
                              tumor_sample = stringr::str_split(string = sample_set$sample[i], pattern = "_vs_", simplify = T)[,1],
                              normal_sample = stringr::str_split(string = sample_set$sample[i], pattern = "_vs_", simplify = T)[,2],
                              caller = "varscan",
                              mut_type = mut_type)
  # Get the caller specific alternate read depth, total read depth, and VAF into common format
  varscan_vaf_dt <- get_vaf(vcf_obj = varscan_gr,
                           caller = "varscan",
                           mut_type = mut_type)
  # Add the common metadata to the GR obj
  S4Vectors::mcols(varscan_gr) <- c(S4Vectors::mcols(varscan_gr), varscan_vaf_dt)
  # Add prefix to column names to keep uniqueness
  colnames(mcols(varscan_gr))[8:length(mcols(varscan_gr))] <- stringr::str_c("VARSCAN_",
                                                                             colnames(mcols(varscan_gr))[8:length(mcols(varscan_gr))])
  
  # Edge case: include CaVEMan for SNV consensus
  if(mut_type == "snv" & length(caveman_vcf_input) == 1) {
    caveman_gr <- read_vcf_file(vcf_file_path = paste0(vcf_file_dir, caveman_vcf),
                                tumor_sample = stringr::str_split(string = sample_set$sample[i], pattern = "_vs_", simplify = T)[,1],
                                normal_sample = stringr::str_split(string = sample_set$sample[i], pattern = "_vs_", simplify = T)[,2],
                                caller = "caveman",
                                mut_type = "snv")
    # Get the caller specific alternate read depth, total read depth, and VAF into common format
    caveman_vaf_dt <- get_vaf(vcf_obj = caveman_gr,
                              caller = "caveman",
                              mut_type = "snv")
    # Add the common metadata to the GR obj
    S4Vectors::mcols(caveman_gr) <- c(S4Vectors::mcols(caveman_gr), caveman_vaf_dt)
    # Add prefix to column names to keep uniqueness
    colnames(mcols(caveman_gr))[8:length(mcols(caveman_gr))] <- stringr::str_c("CAVEMAN_",
                                                                             colnames(mcols(caveman_gr))[8:length(mcols(caveman_gr))])
    
    # Edge case: include SvABA for InDel consensus
  } else if(mut_type == "indel") {
    svaba_gr <- read_vcf_file(vcf_file_path = paste0(vcf_file_dir, svaba_vcf),
                              tumor_sample = stringr::str_split(string = sample_set$sample[i], pattern = "_vs_", simplify = T)[,1],
                              normal_sample = stringr::str_split(string = sample_set$sample[i], pattern = "_vs_", simplify = T)[,2],
                              caller = "svaba",
                              mut_type = "indel")
    # Get the caller specific alternate read depth, total read depth, and VAF into common format
    svaba_vaf_dt <- get_vaf(vcf_obj = svaba_gr,
                              caller = "svaba",
                              mut_type = "indel")
    # Add the common metadata to the GR obj
    S4Vectors::mcols(svaba_gr) <- c(S4Vectors::mcols(svaba_gr), svaba_vaf_dt)
    # Add prefix to column names to keep uniqueness
    colnames(mcols(svaba_gr))[8:length(mcols(svaba_gr))] <- stringr::str_c("SVABA_",
                                                                           colnames(mcols(svaba_gr))[8:length(mcols(svaba_gr))])
  }
  
  # Merge all the calls into single, sorted GR obj
  if(mut_type == "snv" & length(caveman_vcf_input) == 0) {
    union_gr <- gUtils::grbind(mutect_gr, strelka_gr, varscan_gr)
    
  } else if(mut_type == "snv" & length(caveman_vcf_input) == 1) {
    union_gr <- gUtils::grbind(caveman_gr, mutect_gr, strelka_gr, varscan_gr)
    
  } else if(mut_type == "indel") {
    union_gr <- gUtils::grbind(mutect_gr, strelka_gr, varscan_gr, svaba_gr)
    
  }
  
  
  # Sort and refactor seqinfo for consistency
  union_gr <- GenomicRanges::sort.GenomicRanges(union_gr)
  union_gr <- gr_refactor_seqs(input_gr = union_gr)
  
  # Loop through all merged calls to find and merge consensus records (i.e. same call from multiple callers)
  # Split each input loop by chromosome
  chrom_iter_list <- as.character(union_gr@seqnames@values)
  
  final_union_consensus_gr <- foreach::foreach(x = 1:length(chrom_iter_list), .combine = grbind, .packages = "gUtils") %dopar% {
    
    union_per_chrom_gr <- union_gr %Q% (seqnames == chrom_iter_list[x])
    
    # Create empty GRanges
    union_consensus_gr <- GenomicRanges::GRanges()
    
    for(j in 1:length(union_per_chrom_gr)) {
      
      # Check each record for any overlapping records by matching chrom, start, end, patient, ref, and alt columns
      test_for_overlap <- gUtils::gr.findoverlaps(query = union_per_chrom_gr[j],
                                                  subject = union_per_chrom_gr,
                                                  scol = colnames(mcols(union_per_chrom_gr)),
                                                  by = c("PATIENT", "REF", "ALT"))
      
      # overlap test GR will be 1,2,3, or 4
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
  message("\nAnnotating mutations by nearest gene using ", gene_gtf_file," ...")
  final_union_consensus_gr$nearest_gene <- genes$gene_name[IRanges::nearest(x = final_union_consensus_gr,
                                                                            subject = gUtils::gr.stripstrand(genes))]
  
  # Slim down object before output (rm query.id, subject.id)
  final_union_consensus_gr <- final_union_consensus_gr[,c(-1,-2)]

  # Convert GRanges to DT
  final_union_consensus_dt <- gUtils::gr2dt(x = final_union_consensus_gr)
  
  # Now loop through all union consensus calls to report VAFs in single column
  message("Collating per caller read depth and VAF metrics for all mutations ...")
  final_union_consensus_vaf_metrics_dt <- foreach(i = 1:length(chrom_iter_list), .combine = rrbind, .packages = "gUtils") %dopar% {

    # Split each input loop by chromosome
    final_union_consensus_per_chrom_dt <- final_union_consensus_dt[final_union_consensus_dt$seqnames == chrom_iter_list[i]]

    vaf_info <- data.table(alt_read_depth_combo = rep(NA, nrow(final_union_consensus_per_chrom_dt)),
                           alt_read_depth_mean = rep(NA, nrow(final_union_consensus_per_chrom_dt)),
                           total_depth_combo = rep(NA, nrow(final_union_consensus_per_chrom_dt)),
                           total_depth_mean = rep(NA, nrow(final_union_consensus_per_chrom_dt)),
                           vaf_combo = rep(NA, nrow(final_union_consensus_per_chrom_dt)),
                           vaf_mean = rep(NA, nrow(final_union_consensus_per_chrom_dt)))

    for(j in 1:nrow(final_union_consensus_per_chrom_dt)) {

      # First, determine which callers found the variant to start the decision tree
      per_chrom_mut_record <- final_union_consensus_per_chrom_dt[j,]
      caller_string <- per_chrom_mut_record[,CALLER]

      # Use the string of callers to find correct VAF combo and mean
      if(caller_string == "caveman,mutect,strelka,varscan") {
        # 4 CALLER SNV CONSENSUS
        alt_read_depth_combo_string <- stringr::str_c(per_chrom_mut_record$CAVEMAN_alt_depth, per_chrom_mut_record$MUTECT_alt_depth, per_chrom_mut_record$STRELKA_alt_depth, per_chrom_mut_record$VARSCAN_alt_depth,
                                                      sep = ",")
        mean_alt_read_depth <- round(mean(c(per_chrom_mut_record$CAVEMAN_alt_depth, per_chrom_mut_record$MUTECT_alt_depth, per_chrom_mut_record$STRELKA_alt_depth, per_chrom_mut_record$VARSCAN_alt_depth)),
                                     digits = 0)
        total_read_depth_combo_string <- stringr::str_c(per_chrom_mut_record$CAVEMAN_total_depth, per_chrom_mut_record$MUTECT_total_depth, per_chrom_mut_record$STRELKA_total_depth, per_chrom_mut_record$VARSCAN_total_depth,
                                                        sep = ",")
        mean_total_read_depth <- round(mean(c(per_chrom_mut_record$CAVEMAN_total_depth, per_chrom_mut_record$MUTECT_total_depth, per_chrom_mut_record$STRELKA_total_depth, per_chrom_mut_record$VARSCAN_total_depth)),
                                       digits = 0)
        vaf_combo_string <- str_c(per_chrom_mut_record$CAVEMAN_vaf, per_chrom_mut_record$MUTECT_vaf, per_chrom_mut_record$STRELKA_vaf, per_chrom_mut_record$VARSCAN_vaf,
                                  sep = ",")
        mean_vaf <- round(mean(c(per_chrom_mut_record$CAVEMAN_vaf, per_chrom_mut_record$MUTECT_vaf, per_chrom_mut_record$STRELKA_vaf, per_chrom_mut_record$VARSCAN_vaf)),
                          digits = 4)

      } else if(caller_string == "mutect,strelka,svaba,varscan") {
        # 4 CALLER CONSENSUS
        alt_read_depth_combo_string <- stringr::str_c(per_chrom_mut_record$MUTECT_alt_depth, per_chrom_mut_record$STRELKA_alt_depth, per_chrom_mut_record$SVABA_alt_depth, per_chrom_mut_record$VARSCAN_alt_depth,
                                                      sep = ",")
        mean_alt_read_depth <- round(mean(c(per_chrom_mut_record$MUTECT_alt_depth, per_chrom_mut_record$STRELKA_alt_depth, per_chrom_mut_record$SVABA_alt_depth, per_chrom_mut_record$VARSCAN_alt_depth)),
                                     digits = 0)
        total_read_depth_combo_string <- stringr::str_c(per_chrom_mut_record$MUTECT_total_depth, per_chrom_mut_record$STRELKA_total_depth,  per_chrom_mut_record$SVABA_total_depth, per_chrom_mut_record$VARSCAN_total_depth,
                                                        sep = ",")
        mean_total_read_depth <- round(mean(c(per_chrom_mut_record$MUTECT_total_depth, per_chrom_mut_record$STRELKA_total_depth, per_chrom_mut_record$SVABA_total_depth, per_chrom_mut_record$VARSCAN_total_depth)),
                                       digits = 0)
        vaf_combo_string <- str_c(per_chrom_mut_record$MUTECT_vaf, per_chrom_mut_record$STRELKA_vaf, per_chrom_mut_record$SVABA_vaf, per_chrom_mut_record$VARSCAN_vaf,
                                  sep = ",")
        mean_vaf <- round(mean(c(per_chrom_mut_record$MUTECT_vaf, per_chrom_mut_record$STRELKA_vaf, per_chrom_mut_record$SVABA_vaf, per_chrom_mut_record$VARSCAN_vaf)),
                          digits = 4)

      } else if(caller_string == "caveman,mutect,strelka") {
        # 3 CALLER CONSENSUS
        alt_read_depth_combo_string <- stringr::str_c(per_chrom_mut_record$CAVEMAN_alt_depth,per_chrom_mut_record$MUTECT_alt_depth, per_chrom_mut_record$STRELKA_alt_depth,
                                                      sep = ",")
        mean_alt_read_depth <- round(mean(c(per_chrom_mut_record$CAVEMAN_alt_depth, per_chrom_mut_record$MUTECT_alt_depth, per_chrom_mut_record$STRELKA_alt_depth)),
                                     digits = 0)
        total_read_depth_combo_string <- stringr::str_c(per_chrom_mut_record$CAVEMAN_total_depth, per_chrom_mut_record$MUTECT_total_depth, per_chrom_mut_record$STRELKA_total_depth,
                                                        sep = ",")
        mean_total_read_depth <- round(mean(c(per_chrom_mut_record$CAVEMAN_total_depth, per_chrom_mut_record$MUTECT_total_depth, per_chrom_mut_record$STRELKA_total_depth)),
                                       digits = 0)
        vaf_combo_string <- str_c(per_chrom_mut_record$CAVEMAN_vaf, per_chrom_mut_record$MUTECT_vaf, per_chrom_mut_record$STRELKA_vaf,
                                  sep = ",")
        mean_vaf <- round(mean(c(per_chrom_mut_record$CAVEMAN_vaf, per_chrom_mut_record$MUTECT_vaf, per_chrom_mut_record$STRELKA_vaf)),
                          digits = 4)

      } else if(caller_string == "caveman,mutect,varscan") {
        # 3 CALLER CONSENSUS
        alt_read_depth_combo_string <- stringr::str_c(per_chrom_mut_record$CAVEMAN_alt_depth, per_chrom_mut_record$MUTECT_alt_depth, per_chrom_mut_record$VARSCAN_alt_depth,
                                                      sep = ",")
        mean_alt_read_depth <- round(mean(c(per_chrom_mut_record$CAVEMAN_alt_depth, per_chrom_mut_record$MUTECT_alt_depth, per_chrom_mut_record$VARSCAN_alt_depth)),
                                     digits = 0)
        total_read_depth_combo_string <- stringr::str_c(per_chrom_mut_record$CAVEMAN_total_depth, per_chrom_mut_record$MUTECT_total_depth, per_chrom_mut_record$VARSCAN_total_depth,
                                                        sep = ",")
        mean_total_read_depth <- round(mean(c(per_chrom_mut_record$CAVEMAN_total_depth, per_chrom_mut_record$MUTECT_total_depth, per_chrom_mut_record$VARSCAN_total_depth)),
                                       digits = 0)
        vaf_combo_string <- str_c(per_chrom_mut_record$CAVEMAN_vaf, per_chrom_mut_record$MUTECT_vaf, per_chrom_mut_record$VARSCAN_vaf,
                                  sep = ",")
        mean_vaf <- round(mean(c(per_chrom_mut_record$CAVEMAN_vaf, per_chrom_mut_record$MUTECT_vaf, per_chrom_mut_record$VARSCAN_vaf)),
                          digits = 4)

      } else if(caller_string == "caveman,strelka,varscan") {
        # 3 CALLER CONSENSUS
        alt_read_depth_combo_string <- stringr::str_c(per_chrom_mut_record$CAVEMAN_alt_depth, per_chrom_mut_record$STRELKA_alt_depth, per_chrom_mut_record$VARSCAN_alt_depth,
                                                      sep = ",")
        mean_alt_read_depth <- round(mean(c(per_chrom_mut_record$CAVEMAN_alt_depth, per_chrom_mut_record$STRELKA_alt_depth, per_chrom_mut_record$VARSCAN_alt_depth)),
                                     digits = 0)
        total_read_depth_combo_string <- stringr::str_c(per_chrom_mut_record$CAVEMAN_total_depth, per_chrom_mut_record$STRELKA_total_depth, per_chrom_mut_record$VARSCAN_total_depth,
                                                        sep = ",")
        mean_total_read_depth <- round(mean(c(per_chrom_mut_record$CAVEMAN_total_depth, per_chrom_mut_record$STRELKA_total_depth, per_chrom_mut_record$VARSCAN_total_depth)),
                                       digits = 0)
        vaf_combo_string <- str_c(per_chrom_mut_record$CAVEMAN_vaf, per_chrom_mut_record$STRELKA_vaf, per_chrom_mut_record$VARSCAN_vaf,
                                  sep = ",")
        mean_vaf <- round(mean(c(per_chrom_mut_record$CAVEMAN_vaf, per_chrom_mut_record$STRELKA_vaf, per_chrom_mut_record$VARSCAN_vaf)),
                          digits = 4)

      } else if(caller_string == "mutect,strelka,svaba") {
        # 3 CALLER CONSENSUS
        alt_read_depth_combo_string <- stringr::str_c(per_chrom_mut_record$MUTECT_alt_depth, per_chrom_mut_record$STRELKA_alt_depth, per_chrom_mut_record$SVABA_alt_depth,
                                                      sep = ",")
        mean_alt_read_depth <- round(mean(c(per_chrom_mut_record$MUTECT_alt_depth, per_chrom_mut_record$STRELKA_alt_depth, per_chrom_mut_record$SVABA_alt_depth)),
                                     digits = 0)
        total_read_depth_combo_string <- stringr::str_c(per_chrom_mut_record$MUTECT_total_depth, per_chrom_mut_record$STRELKA_total_depth, per_chrom_mut_record$SVABA_total_depth,
                                                        sep = ",")
        mean_total_read_depth <- round(mean(c(per_chrom_mut_record$MUTECT_total_depth, per_chrom_mut_record$STRELKA_total_depth, per_chrom_mut_record$SVABA_total_depth)),
                                       digits = 0)
        vaf_combo_string <- str_c(per_chrom_mut_record$MUTECT_vaf, per_chrom_mut_record$STRELKA_vaf, per_chrom_mut_record$SVABA_vaf,
                                  sep = ",")
        mean_vaf <- round(mean(c(per_chrom_mut_record$MUTECT_vaf, per_chrom_mut_record$STRELKA_vaf, per_chrom_mut_record$SVABA_vaf)),
                          digits = 4)

      } else if(caller_string == "mutect,strelka,varscan") {
        # 3 CALLER CONSENSUS
        alt_read_depth_combo_string <- stringr::str_c(per_chrom_mut_record$MUTECT_alt_depth, per_chrom_mut_record$STRELKA_alt_depth, per_chrom_mut_record$VARSCAN_alt_depth,
                                                      sep = ",")
        mean_alt_read_depth <- round(mean(c(per_chrom_mut_record$MUTECT_alt_depth, per_chrom_mut_record$STRELKA_alt_depth, per_chrom_mut_record$VARSCAN_alt_depth)),
                                     digits = 0)
        total_read_depth_combo_string <- stringr::str_c(per_chrom_mut_record$MUTECT_total_depth, per_chrom_mut_record$STRELKA_total_depth, per_chrom_mut_record$VARSCAN_total_depth,
                                                        sep = ",")
        mean_total_read_depth <- round(mean(c(per_chrom_mut_record$MUTECT_total_depth, per_chrom_mut_record$STRELKA_total_depth, per_chrom_mut_record$VARSCAN_total_depth)),
                                       digits = 0)
        vaf_combo_string <- str_c(per_chrom_mut_record$MUTECT_vaf, per_chrom_mut_record$STRELKA_vaf, per_chrom_mut_record$VARSCAN_vaf,
                                  sep = ",")
        mean_vaf <- round(mean(c(per_chrom_mut_record$MUTECT_vaf, per_chrom_mut_record$STRELKA_vaf, per_chrom_mut_record$VARSCAN_vaf)),
                          digits = 4)

      } else if(caller_string == "mutect,svaba,varscan") {
        # 3 CALLER CONSENSUS
        alt_read_depth_combo_string <- stringr::str_c(per_chrom_mut_record$MUTECT_alt_depth, per_chrom_mut_record$SVABA_alt_depth, per_chrom_mut_record$VARSCAN_alt_depth,
                                                      sep = ",")
        mean_alt_read_depth <- round(mean(c(per_chrom_mut_record$MUTECT_alt_depth, per_chrom_mut_record$SVABA_alt_depth, per_chrom_mut_record$VARSCAN_alt_depth)),
                                     digits = 0)
        total_read_depth_combo_string <- stringr::str_c(per_chrom_mut_record$MUTECT_total_depth, per_chrom_mut_record$SVABA_total_depth, per_chrom_mut_record$VARSCAN_total_depth,
                                                        sep = ",")
        mean_total_read_depth <- round(mean(c(per_chrom_mut_record$MUTECT_total_depth, per_chrom_mut_record$SVABA_total_depth, per_chrom_mut_record$VARSCAN_total_depth)),
                                       digits = 0)
        vaf_combo_string <- str_c(per_chrom_mut_record$MUTECT_vaf, per_chrom_mut_record$SVABA_vaf, per_chrom_mut_record$VARSCAN_vaf,
                                  sep = ",")
        mean_vaf <- round(mean(c(per_chrom_mut_record$MUTECT_vaf, per_chrom_mut_record$SVABA_vaf, per_chrom_mut_record$VARSCAN_vaf)),
                          digits = 4)

      } else if(caller_string == "strelka,svaba,varscan") {
        # 3 CALLER CONSENSUS
        alt_read_depth_combo_string <- stringr::str_c(per_chrom_mut_record$STRELKA_alt_depth, per_chrom_mut_record$SVABA_alt_depth, per_chrom_mut_record$VARSCAN_alt_depth,
                                                      sep = ",")
        mean_alt_read_depth <- round(mean(c(per_chrom_mut_record$STRELKA_alt_depth, per_chrom_mut_record$SVABA_alt_depth, per_chrom_mut_record$VARSCAN_alt_depth)),
                                     digits = 0)
        total_read_depth_combo_string <- stringr::str_c(per_chrom_mut_record$STRELKA_total_depth, per_chrom_mut_record$SVABA_total_depth, per_chrom_mut_record$VARSCAN_total_depth,
                                                        sep = ",")
        mean_total_read_depth <- round(mean(c(per_chrom_mut_record$STRELKA_total_depth, per_chrom_mut_record$SVABA_total_depth, per_chrom_mut_record$VARSCAN_total_depth)),
                                       digits = 0)
        vaf_combo_string <- str_c(per_chrom_mut_record$STRELKA_vaf, per_chrom_mut_record$SVABA_vaf, per_chrom_mut_record$VARSCAN_vaf,
                                  sep = ",")
        mean_vaf <- round(mean(c(per_chrom_mut_record$STRELKA_vaf, per_chrom_mut_record$SVABA_vaf, per_chrom_mut_record$VARSCAN_vaf)),
                          digits = 4)

      } else if(caller_string == "caveman,mutect") {
        # 3 CALLER CONSENSUS
        alt_read_depth_combo_string <- stringr::str_c(per_chrom_mut_record$CAVEMAN_alt_depth, per_chrom_mut_record$MUTECT_alt_depth,
                                                      sep = ",")
        mean_alt_read_depth <- round(mean(c(per_chrom_mut_record$CAVEMAN_alt_depth, per_chrom_mut_record$MUTECT_alt_depth)),
                                     digits = 0)
        total_read_depth_combo_string <- stringr::str_c(per_chrom_mut_record$CAVEMAN_total_depth, per_chrom_mut_record$MUTECT_total_depth,
                                                        sep = ",")
        mean_total_read_depth <- round(mean(c(per_chrom_mut_record$CAVEMAN_total_depth, per_chrom_mut_record$MUTECT_total_depth)),
                                       digits = 0)
        vaf_combo_string <- str_c(per_chrom_mut_record$CAVEMAN_vaf, per_chrom_mut_record$MUTECT_vaf,
                                  sep = ",")
        mean_vaf <- round(mean(c(per_chrom_mut_record$CAVEMAN_vaf, per_chrom_mut_record$MUTECT_vaf)),
                          digits = 4)

      } else if(caller_string == "caveman,strelka") {
        # 3 CALLER CONSENSUS
        alt_read_depth_combo_string <- stringr::str_c(per_chrom_mut_record$CAVEMAN_alt_depth, per_chrom_mut_record$STRELKA_alt_depth,
                                                      sep = ",")
        mean_alt_read_depth <- round(mean(c(per_chrom_mut_record$CAVEMAN_alt_depth, per_chrom_mut_record$STRELKA_alt_depth)),
                                     digits = 0)
        total_read_depth_combo_string <- stringr::str_c(per_chrom_mut_record$CAVEMAN_total_depth, per_chrom_mut_record$STRELKA_total_depth,
                                                        sep = ",")
        mean_total_read_depth <- round(mean(c(per_chrom_mut_record$CAVEMAN_total_depth, per_chrom_mut_record$STRELKA_total_depth)),
                                       digits = 0)
        vaf_combo_string <- str_c(per_chrom_mut_record$CAVEMAN_vaf, per_chrom_mut_record$STRELKA_vaf,
                                  sep = ",")
        mean_vaf <- round(mean(c(per_chrom_mut_record$CAVEMAN_vaf, per_chrom_mut_record$STRELKA_vaf)),
                          digits = 4)

      } else if(caller_string == "caveman,varscan") {
        # 3 CALLER CONSENSUS
        alt_read_depth_combo_string <- stringr::str_c(per_chrom_mut_record$CAVEMAN_alt_depth, per_chrom_mut_record$VARSCAN_alt_depth,
                                                      sep = ",")
        mean_alt_read_depth <- round(mean(c(per_chrom_mut_record$CAVEMAN_alt_depth, per_chrom_mut_record$VARSCAN_alt_depth)),
                                     digits = 0)
        total_read_depth_combo_string <- stringr::str_c(per_chrom_mut_record$CAVEMAN_total_depth, per_chrom_mut_record$VARSCAN_total_depth,
                                                        sep = ",")
        mean_total_read_depth <- round(mean(c(per_chrom_mut_record$CAVEMAN_total_depth, per_chrom_mut_record$VARSCAN_total_depth)),
                                       digits = 0)
        vaf_combo_string <- str_c(per_chrom_mut_record$CAVEMAN_vaf, per_chrom_mut_record$VARSCAN_vaf,
                                  sep = ",")
        mean_vaf <- round(mean(c(per_chrom_mut_record$CAVEMAN_vaf, per_chrom_mut_record$VARSCAN_vaf)),
                          digits = 4)

      } else if(caller_string == "mutect,strelka") {
        # 2 CALLER CONSENSUS
        alt_read_depth_combo_string <- stringr::str_c(per_chrom_mut_record$MUTECT_alt_depth, per_chrom_mut_record$STRELKA_alt_depth,
                                                      sep = ",")
        mean_alt_read_depth <- round(mean(c(per_chrom_mut_record$MUTECT_alt_depth, per_chrom_mut_record$STRELKA_alt_depth)),
                                     digits = 0)
        total_read_depth_combo_string <- stringr::str_c(per_chrom_mut_record$MUTECT_total_depth, per_chrom_mut_record$STRELKA_total_depth,
                                                        sep = ",")
        mean_total_read_depth <- round(mean(c(per_chrom_mut_record$MUTECT_total_depth, per_chrom_mut_record$STRELKA_total_depth)),
                                       digits = 0)
        vaf_combo_string <- str_c(per_chrom_mut_record$MUTECT_vaf, per_chrom_mut_record$STRELKA_vaf,
                                  sep = ",")
        mean_vaf <- round(mean(c(per_chrom_mut_record$MUTECT_vaf, per_chrom_mut_record$STRELKA_vaf)),
                          digits = 4)

      } else if(caller_string == "mutect,svaba") {
        # 2 CALLER CONSENSUS
        alt_read_depth_combo_string <- stringr::str_c(per_chrom_mut_record$MUTECT_alt_depth, per_chrom_mut_record$SVABA_alt_depth,
                                                      sep = ",")
        mean_alt_read_depth <- round(mean(c(per_chrom_mut_record$MUTECT_alt_depth, per_chrom_mut_record$SVABA_alt_depth)),
                                     digits = 0)
        total_read_depth_combo_string <- stringr::str_c(per_chrom_mut_record$MUTECT_total_depth, per_chrom_mut_record$SVABA_total_depth,
                                                        sep = ",")
        mean_total_read_depth <- round(mean(c(per_chrom_mut_record$MUTECT_total_depth, per_chrom_mut_record$SVABA_total_depth)),
                                       digits = 0)
        vaf_combo_string <- str_c(per_chrom_mut_record$MUTECT_vaf, per_chrom_mut_record$SVABA_vaf,
                                  sep = ",")
        mean_vaf <- round(mean(c(per_chrom_mut_record$MUTECT_vaf, per_chrom_mut_record$SVABA_vaf)),
                          digits = 4)

      } else if(caller_string == "mutect,varscan") {
        # 2 CALLER CONSENSUS
        alt_read_depth_combo_string <- stringr::str_c(per_chrom_mut_record$MUTECT_alt_depth, per_chrom_mut_record$VARSCAN_alt_depth,
                                                      sep = ",")
        mean_alt_read_depth <- round(mean(c(per_chrom_mut_record$MUTECT_alt_depth, per_chrom_mut_record$VARSCAN_alt_depth)),
                                     digits = 0)
        total_read_depth_combo_string <- stringr::str_c(per_chrom_mut_record$MUTECT_total_depth, per_chrom_mut_record$VARSCAN_total_depth,
                                                        sep = ",")
        mean_total_read_depth <- round(mean(c(per_chrom_mut_record$MUTECT_total_depth, per_chrom_mut_record$VARSCAN_total_depth)),
                                       digits = 0)
        vaf_combo_string <- str_c(per_chrom_mut_record$MUTECT_vaf, per_chrom_mut_record$VARSCAN_vaf,
                                  sep = ",")
        mean_vaf <- round(mean(c(per_chrom_mut_record$MUTECT_vaf, per_chrom_mut_record$VARSCAN_vaf)),
                          digits = 4)

      } else if(caller_string == "strelka,svaba") {
        # 2 CALLER CONSENSUS
        alt_read_depth_combo_string <- stringr::str_c(per_chrom_mut_record$STRELKA_alt_depth, per_chrom_mut_record$SVABA_alt_depth,
                                                      sep = ",")
        mean_alt_read_depth <- round(mean(c(per_chrom_mut_record$STRELKA_alt_depth, per_chrom_mut_record$SVABA_alt_depth)),
                                     digits = 0)
        total_read_depth_combo_string <- stringr::str_c(per_chrom_mut_record$STRELKA_total_depth, per_chrom_mut_record$SVABA_total_depth,
                                                        sep = ",")
        mean_total_read_depth <- round(mean(c(per_chrom_mut_record$STRELKA_total_depth, per_chrom_mut_record$SVABA_total_depth)),
                                       digits = 0)
        vaf_combo_string <- str_c(per_chrom_mut_record$STRELKA_vaf, per_chrom_mut_record$SVABA_vaf,
                                  sep = ",")
        mean_vaf <- round(mean(c(per_chrom_mut_record$STRELKA_vaf, per_chrom_mut_record$SVABA_vaf)),
                          digits = 4)

      } else if(caller_string == "strelka,varscan") {
        # 2 CALLER CONSENSUS
        alt_read_depth_combo_string <- stringr::str_c(per_chrom_mut_record$STRELKA_alt_depth, per_chrom_mut_record$VARSCAN_alt_depth,
                                                      sep = ",")
        mean_alt_read_depth <- round(mean(c(per_chrom_mut_record$STRELKA_alt_depth, per_chrom_mut_record$VARSCAN_alt_depth)),
                                     digits = 0)
        total_read_depth_combo_string <- stringr::str_c(per_chrom_mut_record$STRELKA_total_depth, per_chrom_mut_record$VARSCAN_total_depth,
                                                        sep = ",")
        mean_total_read_depth <- round(mean(c(per_chrom_mut_record$STRELKA_total_depth, per_chrom_mut_record$VARSCAN_total_depth)),
                                       digits = 0)
        vaf_combo_string <- str_c(per_chrom_mut_record$STRELKA_vaf, per_chrom_mut_record$VARSCAN_vaf,
                                  sep = ",")
        mean_vaf <- round(mean(c(per_chrom_mut_record$STRELKA_vaf, per_chrom_mut_record$VARSCAN_vaf)),
                          digits = 4)

      } else if(caller_string == "svaba,varscan") {
        # 2 CALLER CONSENSUS
        alt_read_depth_combo_string <- stringr::str_c(per_chrom_mut_record$SVABA_alt_depth, per_chrom_mut_record$VARSCAN_alt_depth,
                                                      sep = ",")
        mean_alt_read_depth <- round(mean(c(per_chrom_mut_record$SVABA_alt_depth, per_chrom_mut_record$VARSCAN_alt_depth)),
                                     digits = 0)
        total_read_depth_combo_string <- stringr::str_c(per_chrom_mut_record$SVABA_total_depth, per_chrom_mut_record$VARSCAN_total_depth,
                                                        sep = ",")
        mean_total_read_depth <- round(mean(c(per_chrom_mut_record$SVABA_total_depth, per_chrom_mut_record$VARSCAN_total_depth)),
                                       digits = 0)
        vaf_combo_string <- str_c(per_chrom_mut_record$SVABA_vaf, per_chrom_mut_record$VARSCAN_vaf,
                                  sep = ",")
        mean_vaf <- round(mean(c(per_chrom_mut_record$SVABA_vaf, per_chrom_mut_record$VARSCAN_vaf)),
                          digits = 4)

      } else if(caller_string == "caveman") {
        # 3 CALLER CONSENSUS
        alt_read_depth_combo_string <- per_chrom_mut_record$CAVEMAN_alt_depth
        mean_alt_read_depth <- per_chrom_mut_record$CAVEMAN_alt_depth
        total_read_depth_combo_string <- per_chrom_mut_record$CAVEMAN_total_depth
        mean_total_read_depth <- per_chrom_mut_record$CAVEMAN_total_depth
        vaf_combo_string <- per_chrom_mut_record$CAVEMAN_vaf
        mean_vaf <- per_chrom_mut_record$CAVEMAN_vaf

      } else if(caller_string == "mutect") {
        # 1 CALLER SINGLETON
        alt_read_depth_combo_string <- per_chrom_mut_record$MUTECT_alt_depth
        mean_alt_read_depth <- per_chrom_mut_record$MUTECT_alt_depth
        total_read_depth_combo_string <- per_chrom_mut_record$MUTECT_total_depth
        mean_total_read_depth <- per_chrom_mut_record$MUTECT_total_depth
        vaf_combo_string <- per_chrom_mut_record$MUTECT_vaf
        mean_vaf <- per_chrom_mut_record$MUTECT_vaf

      } else if(caller_string == "strelka") {
        # 1 CALLER SINGLETON
        alt_read_depth_combo_string <- per_chrom_mut_record$STRELKA_alt_depth
        mean_alt_read_depth <- per_chrom_mut_record$STRELKA_alt_depth
        total_read_depth_combo_string <- per_chrom_mut_record$STRELKA_total_depth
        mean_total_read_depth <- per_chrom_mut_record$STRELKA_total_depth
        vaf_combo_string <- per_chrom_mut_record$STRELKA_vaf
        mean_vaf <- per_chrom_mut_record$STRELKA_vaf

      } else if(caller_string == "svaba") {
        # 1 CALLER SINGLETON
        alt_read_depth_combo_string <- per_chrom_mut_record$SVABA_alt_depth
        mean_alt_read_depth <- per_chrom_mut_record$SVABA_alt_depth
        total_read_depth_combo_string <- per_chrom_mut_record$SVABA_total_depth
        mean_total_read_depth <- per_chrom_mut_record$SVABA_total_depth
        vaf_combo_string <- per_chrom_mut_record$SVABA_vaf
        mean_vaf <- per_chrom_mut_record$SVABA_vaf

      } else if(caller_string == "varscan") {
        # 1 CALLER SINGLETON
        alt_read_depth_combo_string <- per_chrom_mut_record$VARSCAN_alt_depth
        mean_alt_read_depth <- per_chrom_mut_record$VARSCAN_alt_depth
        total_read_depth_combo_string <- per_chrom_mut_record$VARSCAN_total_depth
        mean_total_read_depth <- per_chrom_mut_record$VARSCAN_total_depth
        vaf_combo_string <- per_chrom_mut_record$VARSCAN_vaf
        mean_vaf <- per_chrom_mut_record$VARSCAN_vaf
      }

      # Add all VAF info into output DT record
      vaf_info$alt_read_depth_combo[j] <- alt_read_depth_combo_string
      vaf_info$alt_read_depth_mean[j] <- mean_alt_read_depth

      vaf_info$total_depth_combo[j] <- total_read_depth_combo_string
      vaf_info$total_depth_mean[j] <- mean_total_read_depth

      vaf_info$vaf_combo[j] <- vaf_combo_string
      vaf_info$vaf_mean[j] <- mean_vaf
    }

    # print foreach output to conclude parallel run
    vaf_info
  }
  
  # bind the VAF metrics with the main mutation DT
  final_union_consensus_vaf_dt <- cbind(final_union_consensus_dt, final_union_consensus_vaf_metrics_dt)

  # Write the output
  message("Writing output .txt files to ", output_dir, " ...")
  data.table::fwrite(x = final_union_consensus_vaf_dt,
                     file = paste0(output_dir, "/", sample_set$sample[i], ".hq.union.consensus.somatic.", mut_type, ".txt.gz"),
                     quote = FALSE,
                     sep = "\t",
                     row.names = FALSE,
                     col.names = TRUE,
                     logical01 = FALSE,
                     na = "NA",
                     nThread = threads)
  
  message("D O N E ...\n")
}
