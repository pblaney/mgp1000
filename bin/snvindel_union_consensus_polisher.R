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

#' Calculate VAF metrics for tumor and normal samples
#' @param mut_record Single row data.table/data.frame of mutation record with caller-specific columns
#' @param caller_string String of callers separated by commas (e.g., "mutect,strelka,varscan")
#' @param mut_type Type of mutation: "snv" or "indel"
#' @return List with tumor and normal metrics maintaining original column names
calculate_vaf_metrics <- function(mut_record, caller_string, mut_type) {
  
  # Split caller string into individual callers
  callers <- stringr::str_split(caller_string, ",")[[1]]
  
  # Initialize result vectors for TUMOR
  tumor_alt_depths <- c()
  tumor_total_depths <- c()
  tumor_vafs <- c()
  
  # Initialize result vectors for NORMAL
  normal_alt_depths <- c()
  normal_total_depths <- c()
  normal_vafs <- c()
  
  # Process each caller
  for(caller in callers) {
    caller_upper <- stringr::str_to_upper(caller)
    
    # Check if this caller has data for this record
    check_col <- paste0(caller_upper, "_DP_TUMOR")
    if(!check_col %in% colnames(mut_record) || is.na(mut_record[[check_col]])) {
      next
    }
    
    # ========== MUTECT ==========
    if(caller == "mutect") {
      # TUMOR
      tumor_alt_depths <- c(tumor_alt_depths, mut_record[[paste0(caller_upper, "_AD_ALT_TUMOR")]])
      tumor_total_depths <- c(tumor_total_depths, mut_record[[paste0(caller_upper, "_DP_TUMOR")]])
      tumor_vafs <- c(tumor_vafs, mut_record[[paste0(caller_upper, "_AF_TUMOR")]])
      
      # NORMAL
      normal_alt_depths <- c(normal_alt_depths, mut_record[[paste0(caller_upper, "_AD_ALT_NORMAL")]])
      normal_total_depths <- c(normal_total_depths, mut_record[[paste0(caller_upper, "_DP_NORMAL")]])
      normal_vafs <- c(normal_vafs, 
                      mut_record[[paste0(caller_upper, "_AD_ALT_NORMAL")]] / 
                      mut_record[[paste0(caller_upper, "_DP_NORMAL")]])
      
    # ========== VARSCAN ==========
    } else if(caller == "varscan") {
      # TUMOR
      tumor_alt_depths <- c(tumor_alt_depths, mut_record[[paste0(caller_upper, "_AD_TUMOR")]])
      tumor_total_depths <- c(tumor_total_depths, mut_record[[paste0(caller_upper, "_DP_TUMOR")]])
      freq_val_tumor <- as.numeric(stringr::str_remove(mut_record[[paste0(caller_upper, "_FREQ_TUMOR")]], "%")) / 100
      tumor_vafs <- c(tumor_vafs, freq_val_tumor)
      
      # NORMAL
      normal_alt_depths <- c(normal_alt_depths, mut_record[[paste0(caller_upper, "_AD_NORMAL")]])
      normal_total_depths <- c(normal_total_depths, mut_record[[paste0(caller_upper, "_DP_NORMAL")]])
      freq_val_normal <- as.numeric(stringr::str_remove(mut_record[[paste0(caller_upper, "_FREQ_NORMAL")]], "%")) / 100
      normal_vafs <- c(normal_vafs, freq_val_normal)
      
    # ========== STRELKA SNV ==========
    } else if(caller == "strelka" && mut_type == "snv") {
      alt_allele <- mut_record$ALT
      
      # TUMOR
      alt_col_tumor <- paste0(caller_upper, "_", alt_allele, "U_TIER1_TUMOR")
      if(alt_col_tumor %in% colnames(mut_record) && !is.na(mut_record[[alt_col_tumor]])) {
        tumor_alt_depth <- mut_record[[alt_col_tumor]]
        tumor_total_depth <- mut_record[[paste0(caller_upper, "_DP_TUMOR")]]
        
        tumor_alt_depths <- c(tumor_alt_depths, tumor_alt_depth)
        tumor_total_depths <- c(tumor_total_depths, tumor_total_depth)
        tumor_vafs <- c(tumor_vafs, tumor_alt_depth / tumor_total_depth)
      }
      
      # NORMAL
      alt_col_normal <- paste0(caller_upper, "_", alt_allele, "U_TIER1_NORMAL")
      if(alt_col_normal %in% colnames(mut_record) && !is.na(mut_record[[alt_col_normal]])) {
        normal_alt_depth <- mut_record[[alt_col_normal]]
        normal_total_depth <- mut_record[[paste0(caller_upper, "_DP_NORMAL")]]
        
        normal_alt_depths <- c(normal_alt_depths, normal_alt_depth)
        normal_total_depths <- c(normal_total_depths, normal_total_depth)
        normal_vafs <- c(normal_vafs, normal_alt_depth / normal_total_depth)
      }
      
    # ========== STRELKA INDEL ==========
    } else if(caller == "strelka" && mut_type == "indel") {
      # TUMOR
      tumor_alt_depths <- c(tumor_alt_depths, mut_record[[paste0(caller_upper, "_TIR_TIER1_TUMOR")]])
      tumor_total_depths <- c(tumor_total_depths, mut_record[[paste0(caller_upper, "_DP_TUMOR")]])
      tumor_vafs <- c(tumor_vafs, 
                     mut_record[[paste0(caller_upper, "_TIR_TIER1_TUMOR")]] / 
                     mut_record[[paste0(caller_upper, "_DP_TUMOR")]])
      
      # NORMAL
      normal_alt_depths <- c(normal_alt_depths, mut_record[[paste0(caller_upper, "_TIR_TIER1_NORMAL")]])
      normal_total_depths <- c(normal_total_depths, mut_record[[paste0(caller_upper, "_DP_NORMAL")]])
      normal_vafs <- c(normal_vafs, 
                      mut_record[[paste0(caller_upper, "_TIR_TIER1_NORMAL")]] / 
                      mut_record[[paste0(caller_upper, "_DP_NORMAL")]])
      
    # ========== SVABA INDEL ==========
    } else if(caller == "svaba" && mut_type == "indel") {
      # TUMOR
      tumor_alt_depths <- c(tumor_alt_depths, mut_record[[paste0(caller_upper, "_AD_TUMOR")]])
      tumor_total_depths <- c(tumor_total_depths, mut_record[[paste0(caller_upper, "_DP_TUMOR")]])
      tumor_vafs <- c(tumor_vafs, 
                     mut_record[[paste0(caller_upper, "_AD_TUMOR")]] / 
                     mut_record[[paste0(caller_upper, "_DP_TUMOR")]])
      
      # NORMAL
      normal_alt_depths <- c(normal_alt_depths, mut_record[[paste0(caller_upper, "_AD_NORMAL")]])
      normal_total_depths <- c(normal_total_depths, mut_record[[paste0(caller_upper, "_DP_NORMAL")]])
      normal_vafs <- c(normal_vafs, 
                      mut_record[[paste0(caller_upper, "_AD_NORMAL")]] / 
                      mut_record[[paste0(caller_upper, "_DP_NORMAL")]])
      
    # ========== CAVEMAN SNV ==========
    } else if(caller == "caveman" && mut_type == "snv") {
      alt_allele <- mut_record$ALT
      ref_allele <- mut_record$REF
      
      # TUMOR depth calculations
      fwd_alt_col_t <- paste0(caller_upper, "_F", alt_allele, "Z_TUMOR")
      rev_alt_col_t <- paste0(caller_upper, "_R", alt_allele, "Z_TUMOR")
      fwd_ref_col_t <- paste0(caller_upper, "_F", ref_allele, "Z_TUMOR")
      rev_ref_col_t <- paste0(caller_upper, "_R", ref_allele, "Z_TUMOR")
      
      fwd_alt_t <- if(fwd_alt_col_t %in% colnames(mut_record)) mut_record[[fwd_alt_col_t]] else 0
      rev_alt_t <- if(rev_alt_col_t %in% colnames(mut_record)) mut_record[[rev_alt_col_t]] else 0
      fwd_ref_t <- if(fwd_ref_col_t %in% colnames(mut_record)) mut_record[[fwd_ref_col_t]] else 0
      rev_ref_t <- if(rev_ref_col_t %in% colnames(mut_record)) mut_record[[rev_ref_col_t]] else 0
      
      tumor_alt_depth <- sum(fwd_alt_t, rev_alt_t, na.rm = TRUE)
      tumor_total_depth <- sum(fwd_alt_t, rev_alt_t, fwd_ref_t, rev_ref_t, na.rm = TRUE)
      
      tumor_alt_depths <- c(tumor_alt_depths, tumor_alt_depth)
      tumor_total_depths <- c(tumor_total_depths, tumor_total_depth)
      tumor_vafs <- c(tumor_vafs, mut_record[[paste0(caller_upper, "_PM_TUMOR")]])
      
      # NORMAL depth calculations
      fwd_alt_col_n <- paste0(caller_upper, "_F", alt_allele, "Z_NORMAL")
      rev_alt_col_n <- paste0(caller_upper, "_R", alt_allele, "Z_NORMAL")
      fwd_ref_col_n <- paste0(caller_upper, "_F", ref_allele, "Z_NORMAL")
      rev_ref_col_n <- paste0(caller_upper, "_R", ref_allele, "Z_NORMAL")
      
      fwd_alt_n <- if(fwd_alt_col_n %in% colnames(mut_record)) mut_record[[fwd_alt_col_n]] else 0
      rev_alt_n <- if(rev_alt_col_n %in% colnames(mut_record)) mut_record[[rev_alt_col_n]] else 0
      fwd_ref_n <- if(fwd_ref_col_n %in% colnames(mut_record)) mut_record[[fwd_ref_col_n]] else 0
      rev_ref_n <- if(rev_ref_col_n %in% colnames(mut_record)) mut_record[[rev_ref_col_n]] else 0
      
      normal_alt_depth <- sum(fwd_alt_n, rev_alt_n, na.rm = TRUE)
      normal_total_depth <- sum(fwd_alt_n, rev_alt_n, fwd_ref_n, rev_ref_n, na.rm = TRUE)
      
      normal_alt_depths <- c(normal_alt_depths, normal_alt_depth)
      normal_total_depths <- c(normal_total_depths, normal_total_depth)
      normal_vafs <- c(normal_vafs, normal_alt_depth / normal_total_depth)
    }
  }
  
  # Check if we have valid tumor data
  if(length(tumor_alt_depths) == 0) {
    stop("No valid caller data found in the mutation record")
  }
  
  # Calculate TUMOR consensus metrics (maintaining original column names)
  tumor_metrics <- list(
    alt_read_depth_combo = stringr::str_c(tumor_alt_depths, collapse = ","),
    alt_read_depth_mean = round(mean(tumor_alt_depths, na.rm = TRUE), digits = 0),
    total_depth_combo = stringr::str_c(tumor_total_depths, collapse = ","),
    total_depth_mean = round(mean(tumor_total_depths, na.rm = TRUE), digits = 0),
    vaf_combo = stringr::str_c(tumor_vafs, collapse = ","),
    vaf_mean = round(mean(tumor_vafs, na.rm = TRUE), digits = 4)
  )
  
  # Calculate NORMAL consensus metrics (new columns)
  normal_metrics <- list(
    normal_alt_read_depth_combo = NA_character_,
    normal_alt_read_depth_mean = NA_integer_,
    normal_total_depth_combo = NA_character_,
    normal_total_depth_mean = NA_integer_,
    normal_vaf_combo = NA_character_,
    normal_vaf_mean = NA_real_
  )
  
  if(length(normal_alt_depths) > 0) {
    normal_metrics <- list(
      normal_alt_read_depth_combo = stringr::str_c(normal_alt_depths, collapse = ","),
      normal_alt_read_depth_mean = round(mean(normal_alt_depths, na.rm = TRUE), digits = 0),
      normal_total_depth_combo = stringr::str_c(normal_total_depths, collapse = ","),
      normal_total_depth_mean = round(mean(normal_total_depths, na.rm = TRUE), digits = 0),
      normal_vaf_combo = stringr::str_c(normal_vafs, collapse = ","),
      normal_vaf_mean = round(mean(normal_vafs, na.rm = TRUE), digits = 4)
    )
  }
  
  # Combine and return
  return(c(tumor_metrics, normal_metrics))
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
message("Collecting possible per sample VCFs ...")
mutect_samples <- stringr::str_remove(string = mutect_vcf_input, pattern = ".mutect.somatic.*.vcf.gz") %>%
  tibble::as_tibble_col(column_name = "sample")
message(nrow(mutect_samples), " Mutect VCF(s) detected ...")

strelka_samples <- stringr::str_remove(string = strelka_vcf_input, pattern = ".strelka.somatic.*.vcf.gz") %>%
  tibble::as_tibble_col(column_name = "sample")
message(nrow(strelka_samples), " Strelka VCF(s) detected ...")

varscan_samples <- stringr::str_remove(string = varscan_vcf_input, pattern = ".varscan.somatic.*.vcf.gz") %>%
  tibble::as_tibble_col(column_name = "sample")
message(nrow(varscan_samples), " Varscan VCF(s) detected ...")

# Build sample set using full joins to keep all samples regardless of caller availability
sample_set <- dplyr::full_join(x = mutect_samples, y = strelka_samples, by = "sample") %>%
  dplyr::full_join(y = varscan_samples, by = "sample")

# Add optional callers based on mutation type
if(mut_type == "snv" & length(caveman_vcf_input) > 0) {
  caveman_samples <- stringr::str_remove(string = caveman_vcf_input, pattern = ".caveman.somatic.snv.vcf.gz") %>%
    tibble::as_tibble_col(column_name = "sample")
  message(nrow(caveman_samples), " CaVEMan VCF(s) detected ...")
  sample_set <- dplyr::full_join(sample_set, caveman_samples, by = "sample")
}

if(mut_type == "indel" & length(svaba_vcf_input) > 0) {
  svaba_samples <- stringr::str_remove(string = svaba_vcf_input, pattern = ".svaba.somatic.indel.vcf.gz") %>%
    tibble::as_tibble_col(column_name = "sample")
  message(nrow(svaba_samples), " SvABA VCF(s) detected ...")
  sample_set <- dplyr::full_join(sample_set, svaba_samples, by = "sample")
}

message("\n", nrow(sample_set), " unique sample(s) detected across all callers ...")

# Read in genes
genes <- get_genes_shortcut(gtf_file_path = gene_gtf_file)

# Loop through all per sample triplet VCFs to create consensus mutation set
for(i in 1:nrow(sample_set)) {
  
  message("\n\n----------------------------------------")
  message("Processing sample: ", sample_set$sample[i])
  cat("\n")
  
  # List to store GRanges objects for available callers
  gr_list <- list()
  
  # Check and process Mutect
  mutect_vcf <- mutect_vcf_input[stringr::str_detect(string = mutect_vcf_input, pattern = sample_set$sample[i])]
  if(length(mutect_vcf) > 0) {
    message("  Loading Mutect VCF ...")
    mutect_gr <- read_vcf_file(vcf_file_path = paste0(vcf_file_dir, mutect_vcf),
                               tumor_sample = stringr::str_split(string = sample_set$sample[i], pattern = "_vs_", simplify = T)[,1],
                               normal_sample = stringr::str_split(string = sample_set$sample[i], pattern = "_vs_", simplify = T)[,2],
                               caller = "mutect",
                               mut_type = mut_type)
    mutect_vaf_dt <- get_vaf(vcf_obj = mutect_gr, caller = "mutect", mut_type = mut_type)
    S4Vectors::mcols(mutect_gr) <- c(S4Vectors::mcols(mutect_gr), mutect_vaf_dt)
    colnames(mcols(mutect_gr))[8:length(mcols(mutect_gr))] <- stringr::str_c("MUTECT_",
                                                                             colnames(mcols(mutect_gr))[8:length(mcols(mutect_gr))])
    gr_list[["mutect"]] <- mutect_gr
  }
  
  # Check and process Strelka
  strelka_vcf <- strelka_vcf_input[stringr::str_detect(string = strelka_vcf_input, pattern = sample_set$sample[i])]
  if(length(strelka_vcf) > 0) {
    message("  Loading Strelka VCF ...")
    strelka_gr <- read_vcf_file(vcf_file_path = paste0(vcf_file_dir, strelka_vcf),
                                tumor_sample = stringr::str_split(string = sample_set$sample[i], pattern = "_vs_", simplify = T)[,1],
                                normal_sample = stringr::str_split(string = sample_set$sample[i], pattern = "_vs_", simplify = T)[,2],
                                caller = "strelka",
                                mut_type = mut_type)
    strelka_vaf_dt <- get_vaf(vcf_obj = strelka_gr, caller = "strelka", mut_type = mut_type)
    S4Vectors::mcols(strelka_gr) <- c(S4Vectors::mcols(strelka_gr), strelka_vaf_dt)
    colnames(mcols(strelka_gr))[8:length(mcols(strelka_gr))] <- stringr::str_c("STRELKA_",
                                                                               colnames(mcols(strelka_gr))[8:length(mcols(strelka_gr))])
    gr_list[["strelka"]] <- strelka_gr
  }
  
  # Check and process Varscan
  varscan_vcf <- varscan_vcf_input[stringr::str_detect(string = varscan_vcf_input, pattern = sample_set$sample[i])]
  if(length(varscan_vcf) > 0) {
    message("  Loading Varscan VCF ...")
    varscan_gr <- read_vcf_file(vcf_file_path = paste0(vcf_file_dir, varscan_vcf),
                                tumor_sample = stringr::str_split(string = sample_set$sample[i], pattern = "_vs_", simplify = T)[,1],
                                normal_sample = stringr::str_split(string = sample_set$sample[i], pattern = "_vs_", simplify = T)[,2],
                                caller = "varscan",
                                mut_type = mut_type)
    varscan_vaf_dt <- get_vaf(vcf_obj = varscan_gr, caller = "varscan", mut_type = mut_type)
    S4Vectors::mcols(varscan_gr) <- c(S4Vectors::mcols(varscan_gr), varscan_vaf_dt)
    colnames(mcols(varscan_gr))[8:length(mcols(varscan_gr))] <- stringr::str_c("VARSCAN_",
                                                                               colnames(mcols(varscan_gr))[8:length(mcols(varscan_gr))])
    gr_list[["varscan"]] <- varscan_gr
  }
  
  # Check and process CaVEMan (SNV only)
  if(mut_type == "snv") {
    caveman_vcf <- caveman_vcf_input[stringr::str_detect(string = caveman_vcf_input, pattern = sample_set$sample[i])]
    if(length(caveman_vcf) > 0) {
      message("  Loading CaVEMan VCF ...")
      caveman_gr <- read_vcf_file(vcf_file_path = paste0(vcf_file_dir, caveman_vcf),
                                  tumor_sample = stringr::str_split(string = sample_set$sample[i], pattern = "_vs_", simplify = T)[,1],
                                  normal_sample = stringr::str_split(string = sample_set$sample[i], pattern = "_vs_", simplify = T)[,2],
                                  caller = "caveman",
                                  mut_type = "snv")
      caveman_vaf_dt <- get_vaf(vcf_obj = caveman_gr, caller = "caveman", mut_type = "snv")
      S4Vectors::mcols(caveman_gr) <- c(S4Vectors::mcols(caveman_gr), caveman_vaf_dt)
      colnames(mcols(caveman_gr))[8:length(mcols(caveman_gr))] <- stringr::str_c("CAVEMAN_",
                                                                                 colnames(mcols(caveman_gr))[8:length(mcols(caveman_gr))])
      gr_list[["caveman"]] <- caveman_gr
    }
  }
  
  # Check and process SvABA (InDel only)
  if(mut_type == "indel") {
    svaba_vcf <- svaba_vcf_input[stringr::str_detect(string = svaba_vcf_input, pattern = sample_set$sample[i])]
    if(length(svaba_vcf) > 0) {
      message("  Loading SvABA VCF ...")
      svaba_gr <- read_vcf_file(vcf_file_path = paste0(vcf_file_dir, svaba_vcf),
                                tumor_sample = stringr::str_split(string = sample_set$sample[i], pattern = "_vs_", simplify = T)[,1],
                                normal_sample = stringr::str_split(string = sample_set$sample[i], pattern = "_vs_", simplify = T)[,2],
                                caller = "svaba",
                                mut_type = "indel")
      svaba_vaf_dt <- get_vaf(vcf_obj = svaba_gr, caller = "svaba", mut_type = "indel")
      S4Vectors::mcols(svaba_gr) <- c(S4Vectors::mcols(svaba_gr), svaba_vaf_dt)
      colnames(mcols(svaba_gr))[8:length(mcols(svaba_gr))] <- stringr::str_c("SVABA_",
                                                                             colnames(mcols(svaba_gr))[8:length(mcols(svaba_gr))])
      gr_list[["svaba"]] <- svaba_gr
    }
  }
  
  # Check if we have at least one caller
  if(length(gr_list) == 0) {
    message("  WARNING: No VCF files found for sample ", sample_set$sample[i], ". Skipping ...")
    next
  }
  
  message("  Found ", length(gr_list), " caller(s): ", paste(names(gr_list), collapse = ", "))
  
  # Merge all the calls into single, sorted GR obj
  union_gr <- do.call(gUtils::grbind, gr_list)
  
  # Sort and refactor seqinfo for consistency
  union_gr <- GenomicRanges::sort.GenomicRanges(union_gr)
  union_gr <- gr_refactor_seqs(input_gr = union_gr)
  
  # Loop through all merged calls to find and merge consensus records (i.e. same call from multiple callers)
  # Split each input loop by chromosome
  message("  Finding consensus mutations ...")
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
  message("  Annotating mutations by nearest gene using ", gene_gtf_file," ...")
  final_union_consensus_gr$nearest_gene <- genes$gene_name[IRanges::nearest(x = final_union_consensus_gr,
                                                                            subject = gUtils::gr.stripstrand(genes))]
  
  # Slim down object before output (rm query.id, subject.id)
  final_union_consensus_gr <- final_union_consensus_gr[,c(-1,-2)]
  # Convert GRanges to DT
  final_union_consensus_dt <- gUtils::gr2dt(x = final_union_consensus_gr)
  
  # Now loop through all union consensus calls to report VAFs for both tumor and normal
  message("  Collating per caller read depth and VAF metrics for tumor and normal samples ...")
  final_union_consensus_vaf_metrics_dt <- foreach::foreach(i = 1:length(chrom_iter_list), .combine = rrbind, .packages = "gUtils") %dopar% {
    # Split each input loop by chromosome
    final_union_consensus_per_chrom_dt <- final_union_consensus_dt[final_union_consensus_dt$seqnames == chrom_iter_list[i]]
    
    # Initialize data.table with both tumor and normal columns
    vaf_info <- data.table::data.table(
      alt_read_depth_combo = rep(NA, nrow(final_union_consensus_per_chrom_dt)),
      alt_read_depth_mean = rep(NA, nrow(final_union_consensus_per_chrom_dt)),
      total_depth_combo = rep(NA, nrow(final_union_consensus_per_chrom_dt)),
      total_depth_mean = rep(NA, nrow(final_union_consensus_per_chrom_dt)),
      vaf_combo = rep(NA, nrow(final_union_consensus_per_chrom_dt)),
      vaf_mean = rep(NA, nrow(final_union_consensus_per_chrom_dt)),
      normal_alt_read_depth_combo = rep(NA, nrow(final_union_consensus_per_chrom_dt)),
      normal_alt_read_depth_mean = rep(NA, nrow(final_union_consensus_per_chrom_dt)),
      normal_total_depth_combo = rep(NA, nrow(final_union_consensus_per_chrom_dt)),
      normal_total_depth_mean = rep(NA, nrow(final_union_consensus_per_chrom_dt)),
      normal_vaf_combo = rep(NA, nrow(final_union_consensus_per_chrom_dt)),
      normal_vaf_mean = rep(NA, nrow(final_union_consensus_per_chrom_dt))
    )
    
    for(j in 1:nrow(final_union_consensus_per_chrom_dt)) {
      per_chrom_mut_record <- final_union_consensus_per_chrom_dt[j,]
      caller_string <- per_chrom_mut_record[,CALLER]
      
      metrics <- calculate_vaf_metrics(per_chrom_mut_record, caller_string, mut_type)
      
      # Tumor metrics (original column names)
      vaf_info$alt_read_depth_combo[j] <- metrics$alt_read_depth_combo
      vaf_info$alt_read_depth_mean[j] <- metrics$alt_read_depth_mean
      vaf_info$total_depth_combo[j] <- metrics$total_depth_combo
      vaf_info$total_depth_mean[j] <- metrics$total_depth_mean
      vaf_info$vaf_combo[j] <- metrics$vaf_combo
      vaf_info$vaf_mean[j] <- metrics$vaf_mean
      
      # Normal metrics (new columns)
      vaf_info$normal_alt_read_depth_combo[j] <- metrics$normal_alt_read_depth_combo
      vaf_info$normal_alt_read_depth_mean[j] <- metrics$normal_alt_read_depth_mean
      vaf_info$normal_total_depth_combo[j] <- metrics$normal_total_depth_combo
      vaf_info$normal_total_depth_mean[j] <- metrics$normal_total_depth_mean
      vaf_info$normal_vaf_combo[j] <- metrics$normal_vaf_combo
      vaf_info$normal_vaf_mean[j] <- metrics$normal_vaf_mean
    }
    
    # print foreach output to conclude parallel run
    vaf_info
  }
  
  # bind the VAF metrics with the main mutation DT
  final_union_consensus_vaf_dt <- cbind(final_union_consensus_dt, final_union_consensus_vaf_metrics_dt)
  
  # Write the output
  message("  Writing output .txt files to ", output_dir, " ...")
  data.table::fwrite(x = final_union_consensus_vaf_dt,
                     file = paste0(output_dir, "/", sample_set$sample[i], ".hq.union.consensus.somatic.", mut_type, ".txt.gz"),
                     quote = FALSE,
                     sep = "\t",
                     row.names = FALSE,
                     col.names = TRUE,
                     logical01 = FALSE,
                     na = "NA",
                     nThread = threads)
  
  message("  D O N E for sample ", sample_set$sample[i], "...\n")
}