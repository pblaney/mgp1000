#!/usr/bin/env Rscript
# This script accepts set of SV VCF files as input and generates
# a union consensus SV file
#########################
#####   Libraries   #####
suppressPackageStartupMessages(library(StructuralVariantAnnotation))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))
suppressWarnings(suppressPackageStartupMessages(library(gGnome)))
suppressWarnings(suppressPackageStartupMessages(library(devgru)))
suppressPackageStartupMessages(library(paletteer))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(data.table))
Sys.setenv(DEFAULT_BSGENOME = 'BSgenome.Hsapiens.UCSC.hg38::Hsapiens')
options(scipen = 999)

#########################
#####   Functions   #####

# Original author Jennifer Shelton, https://bitbucket.nygenome.org/projects/WDL/repos/somatic_dna_tools/browse/vcf_to_bedpe.r
# commit ef47c0317e8  14 Jan 2022
## Convert breakpointRanges to BEDPE
vcfToBedpe <- function(vcf) {
  
  sqn <- as.character(seqnames(vcf))
  strand <- as.character(strand(vcf))
  res <- c()
  processed <- c()
  
  for (i in 1:length(vcf)) {
    bnd <- names(vcf)[i]
    partner <- vcf$partner[i]
    partner.idx <- which(names(vcf) == partner)
    
    ## If we don't have exactly one partner, exclude this variant
    if (length(partner.idx) != 1) {
      warning('Missing partner for breakend ', bnd)
      next
    }
    
    ## Check to see if we've already processed this or it's partner
    if (any(c(bnd, partner) %in% processed)) {
      next
    }
    
    ## Combine breakends in single line
    res.i <- c(sqn[i], start(vcf)[i], end(vcf)[i],
               sqn[partner.idx], start(vcf)[partner.idx], end(vcf)[partner.idx],
               'BND', '.', strand[i], strand[partner.idx])
    
    ## Add to result, keep track of processed breakends
    res <- rbind(res, res.i)
    processed <- c(processed, bnd, partner)
  }
  
  ## Add colnames and fill in simple event classifications
  colnames(res) <- c('chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'type', 'score', 'strand1', 'strand2')
  res <- as.data.frame(res, stringsAsFactors = F)
  
  res$type[res$strand1 == '+' & res$strand2 == '-'] <- 'DEL'
  res$type[res$strand1 == '-' & res$strand2 == '+'] <- 'DUP'
  res$type[res$strand1 == '-' & res$strand2 == '-'] <- 'INV'
  res$type[res$strand1 == '+' & res$strand2 == '+'] <- 'INV'
  res$type[res$chr1 != res$chr2] <- 'TRA'
  
  ## Sort by chromosome 
  res <- res[order(factor(res$chr1, levels = levels(seqnames(vcf))), res$start1, res$end1, decreasing = F), ]
  
  ## Simplify coordinates
  res$end1 <- as.numeric(res$start1) + 1
  res$end2 <- as.numeric(res$start2) + 1
  
  colnames(res)[1] <- paste0('#', colnames(res)[1])
  
  return(res)
}

# Process individual caller VCF/TSV to junction object
process_caller <- function(caller_name, input_data, input_type = "vcf", 
                          tum_norm_id, ref_seq_info, temp_files_list) {
  
  message("Processing ", caller_name, " input ...")
  
  # Handle different input types
  if (input_type == "vcf") {
    n_records <- nrow(VariantAnnotation::info(input_data))
  } else if (input_type == "tsv") {
    n_records <- nrow(input_data)
  }
  
  message("  Found ", n_records, " SVs in ", toupper(input_type), " ...")
  
  # Return NA early if no records
  if (n_records == 0) {
    return(list(jnc = NA, temp_files = temp_files_list, has_data = FALSE))
  }
  
  # Caller-specific processing
  if (caller_name == "manta") {
    bp_ranges <- StructuralVariantAnnotation::breakpointRanges(input_data,
                                                               nominalPosition = TRUE,
                                                               inferMissingBreakends = TRUE)
    bedpe <- vcfToBedpe(vcf = bp_ranges)
    
  } else if (caller_name == "svaba") {
    # Adjust SvABA VCF to include END info field
    for (i in 1:n_records) {
      if (input_data@info$SPAN[i] != -1) {
        end_of_sv <- input_data@rowRanges[i]@ranges@start + input_data@info$SPAN[i]
        input_data@info$END[i] <- end_of_sv
      } else {
        input_data@info$END[i] <- NA
      }
    }
    
    bp_ranges <- StructuralVariantAnnotation::breakpointRanges(input_data,
                                                               nominalPosition = TRUE,
                                                               inferMissingBreakends = TRUE)
    bedpe <- vcfToBedpe(vcf = bp_ranges)
    
  } else if (caller_name == "delly") {
    bp_ranges <- StructuralVariantAnnotation::breakpointRanges(input_data,
                                                               nominalPosition = TRUE,
                                                               inferMissingBreakends = TRUE)
    bedpe <- vcfToBedpe(vcf = bp_ranges)
    
  } else if (caller_name == "igcaller") {
    # TSV -> BEDPE
    bedpe <- input_data %>%
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
                   "strand2" = StrandB)
    
    bedpe$type <- dplyr::case_when(
      bedpe$type == "Deletion" ~ "DEL",
      bedpe$type == "Duplication" ~ "DUP",
      bedpe$type == "Gain" ~ "DUP",
      bedpe$type == "Insertion" ~ "INS",
      bedpe$type == "Inversion" ~ "INV",
      bedpe$type == "Translocation" ~ "TRA"
    )
  }
  
  # Check if bedpe is valid and has rows
  if (is.null(bedpe) || nrow(bedpe) == 0) {
    message("  Warning: No valid BEDPE records generated for ", caller_name)
    return(list(jnc = NA, temp_files = temp_files_list, has_data = FALSE))
  }
  
  # Write temp BEDPE
  temp_bedpe_filename <- tempfile(pattern = paste0(tum_norm_id, ".", caller_name),
                                 fileext = ".somatic.sv.temp.bedpe")
  temp_files_list <- append(temp_files_list, temp_bedpe_filename)
  
  write.table(x = bedpe, file = temp_bedpe_filename,
             row.names = FALSE, col.names = TRUE,
             sep = "\t", quote = FALSE)
  
  # BEDPE -> gGnome Junctions
  tryCatch({
    jnc <- gGnome::jJ(rafile = temp_bedpe_filename,
                     chr.convert = FALSE,
                     hg = "hg38",
                     keep.features = TRUE,
                     seqlengths = ref_seq_info)
    
    # Validate junction object
    if (is.null(jnc) || length(jnc$grl) == 0) {
      message("  Warning: No valid junctions created for ", caller_name)
      return(list(jnc = NA, temp_files = temp_files_list, has_data = FALSE))
    }
    
    return(list(jnc = jnc, temp_files = temp_files_list, has_data = TRUE))
    
  }, error = function(e) {
    message("  Error creating junctions for ", caller_name, ": ", e$message)
    return(list(jnc = NA, temp_files = temp_files_list, has_data = FALSE))
  })
}

# Dynamic consensus merging
merge_callers_dynamic <- function(caller_list, caller_names) {
  
  # Build merge arguments dynamically
  merge_args <- list()
  for (i in seq_along(caller_names)) {
    merge_args[[caller_names[i]]] <- caller_list[[i]][, "name"]
  }
  merge_args$pad <- 1000
  
  # Perform merge
  consensus_jnc <- do.call(gGnome::merge, merge_args)
  
  return(consensus_jnc)
}

# Extract caller consensus string dynamically
extract_caller_string <- function(bedpe_pair, caller_names) {
  
  # Safety check
  if (nrow(bedpe_pair) == 0) {
    return("unknown")
  }
  
  caller_presence <- c()
  
  for (caller in caller_names) {
    col_name <- paste0("seen.by.", caller)
    
    # Check if column exists
    if (col_name %in% colnames(bedpe_pair)) {
      # Check if any value is TRUE or > 0
      caller_values <- bedpe_pair[[col_name]]
      
      if (any(!is.na(caller_values) & (caller_values == TRUE | caller_values > 0))) {
        caller_presence <- c(caller_presence, caller)
      }
    }
  }
  
  # Build the caller combo string
  if (length(caller_presence) > 0) {
    return(stringr::str_flatten(caller_presence, collapse = ","))
  } else {
    return("unknown")
  }
}

#########################
#####   Execution   #####

# Accept command line arguments as input
input_args <- commandArgs(trailingOnly = TRUE)

input_file_dir <- input_args[1]
output_dir <- input_args[2]

# Set seqinfo named vector of hg38 chrom names and lengths
ref_seq_info <- c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
                 159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
                 114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
                 58617616, 64444167, 46709983, 50818468, 156040895, 57227415)
names(ref_seq_info) <- paste0("chr", c(seq(1, 22, 1), "X", "Y"))

# Find all input files in directory provided
message("\nScanning path ", input_file_dir, " for input SV files ...")

manta_vcf_input <- list.files(path = input_file_dir,
                              pattern = "*.manta.somatic.sv.final.vcf.gz",
                              full.names = FALSE)

svaba_vcf_input <- list.files(path = input_file_dir,
                              pattern = "*.svaba.somatic.sv.final.vcf.gz",
                              full.names = FALSE)

delly_vcf_input <- list.files(path = input_file_dir,
                              pattern = "*.delly.somatic.sv.final.vcf.gz",
                              full.names = FALSE)

igcaller_tsv_input <- list.files(path = input_file_dir,
                                 pattern = "*.igcaller.oncogenic.rearrangements.tsv",
                                 full.names = FALSE)

# Get sample names and gather input file sets
message("Collecting possible per sample SV files .....")

manta_samples <- stringr::str_remove(string = manta_vcf_input, 
                                     pattern = ".manta.somatic.sv.final.vcf.gz") %>%
  tibble::as_tibble_col(column_name = "sample")
message(nrow(manta_samples), " Manta VCF(s) detected ...")

svaba_samples <- stringr::str_remove(string = svaba_vcf_input, 
                                     pattern = ".svaba.somatic.sv.final.vcf.gz") %>%
  tibble::as_tibble_col(column_name = "sample")
message(nrow(svaba_samples), " SvABA VCF(s) detected ...")

delly_samples <- stringr::str_remove(string = delly_vcf_input, 
                                     pattern = ".delly.somatic.sv.final.vcf.gz") %>%
  tibble::as_tibble_col(column_name = "sample")
message(nrow(delly_samples), " DELLY VCF(s) detected ...")

igcaller_samples <- stringr::str_remove(string = igcaller_tsv_input, 
                                        pattern = ".igcaller.oncogenic.rearrangements.tsv") %>%
  tibble::as_tibble_col(column_name = "sample")
message(nrow(igcaller_samples), " IgCaller TSV(s) detected ...")

# Build sample set using full joins to keep all samples regardless of caller availability
sample_set <- dplyr::full_join(x = manta_samples, y = svaba_samples, by = "sample") %>%
  dplyr::full_join(y = delly_samples, by = "sample") %>%
  dplyr::full_join(y = igcaller_samples, by = "sample")

message("\n", nrow(sample_set), " unique sample(s) detected across all callers ...")

# Loop through all per sample SV files to create consensus SV set
for (i in 1:nrow(sample_set)) {
  
  message("\n\n========================================")
  message("Processing sample: ", sample_set$sample[i])
  cat("\n")
  
  tum_norm_id <- sample_set$sample[i]
  
  # Find paths to input files for this sample
  manta_vcf_path <- list.files(path = input_file_dir,
                               pattern = paste0(sample_set$sample[i], ".manta.somatic.sv.final.vcf.gz"),
                               full.names = TRUE)
  
  svaba_vcf_path <- list.files(path = input_file_dir,
                               pattern = paste0(sample_set$sample[i], ".svaba.somatic.sv.final.vcf.gz"),
                               full.names = TRUE)
  
  delly_vcf_path <- list.files(path = input_file_dir,
                               pattern = paste0(sample_set$sample[i], ".delly.somatic.sv.final.vcf.gz"),
                               full.names = TRUE)
  
  igcaller_tsv_path <- list.files(path = input_file_dir,
                                  pattern = paste0(sample_set$sample[i], ".igcaller.oncogenic.rearrangements.tsv"),
                                  full.names = TRUE)
  
  # Process each caller dynamically
  temp_files <- c()
  caller_jnc_list <- list()
  available_callers <- c()
  
  # Process Manta
  if (length(manta_vcf_path) > 0) {
    manta_vcf <- VariantAnnotation::readVcf(file = manta_vcf_path, genome = "hg38")
    result <- process_caller("manta", manta_vcf, "vcf", tum_norm_id, ref_seq_info, temp_files)
    if (result$has_data) {
      caller_jnc_list[["manta"]] <- result$jnc
      available_callers <- c(available_callers, "manta")
    }
    temp_files <- result$temp_files
  }

  # Process SvABA
  if (length(svaba_vcf_path) > 0) {
    svaba_vcf <- VariantAnnotation::readVcf(file = svaba_vcf_path, genome = "hg38")
    result <- process_caller("svaba", svaba_vcf, "vcf", tum_norm_id, ref_seq_info, temp_files)
    if (result$has_data) {
      caller_jnc_list[["svaba"]] <- result$jnc
      available_callers <- c(available_callers, "svaba")
    }
    temp_files <- result$temp_files
  }

  # Process DELLY
  if (length(delly_vcf_path) > 0) {
    delly_vcf <- VariantAnnotation::readVcf(file = delly_vcf_path, genome = "hg38")
    result <- process_caller("delly", delly_vcf, "vcf", tum_norm_id, ref_seq_info, temp_files)
    if (result$has_data) {
      caller_jnc_list[["delly"]] <- result$jnc
      available_callers <- c(available_callers, "delly")
    }
    temp_files <- result$temp_files
  }

  # Process IgCaller
  # Process IgCaller
  if (length(igcaller_tsv_path) > 0) {
    igcaller_tsv <- readr::read_delim(file = igcaller_tsv_path,
                                      delim = "\t",
                                      col_names = TRUE,
                                      show_col_types = FALSE)
    
    # Store original count for reporting
    original_count <- nrow(igcaller_tsv)
    
    # Apply base filtering (existing filter)
    igcaller_tsv <- igcaller_tsv %>%
      dplyr::filter(Score >= 10 & `Reads in normal` <= 1 & `Count in PoN` <= 1)
    
    base_filtered_count <- nrow(igcaller_tsv)
    
    # Identify sex chromosome translocations
    # A translocation involves different chromosomes, and at least one is chrX or chrY
    igcaller_tsv <- igcaller_tsv %>%
      dplyr::mutate(
        is_sex_chr_translocation = (ChrA != ChrB) &  # Different chromosomes (translocation)
                                  ((ChrA %in% c("chrX", "chrY")) | (ChrB %in% c("chrX", "chrY")))
      )
    
    # Apply stringent filtering for sex chromosome translocations
    # Keep sex chr translocations only if they meet strict criteria
    # Keep all non-sex chr translocations if they pass base filter
    igcaller_tsv <- igcaller_tsv %>%
      dplyr::filter(
        # Either it's NOT a sex chromosome translocation (keep if passed base filter)
        !is_sex_chr_translocation |
        # OR it IS a sex chr translocation but meets stringent criteria
        (is_sex_chr_translocation & Score >= 15 & `Reads in normal` == 0 & `Count in PoN` == 0)
      ) %>%
      dplyr::select(-is_sex_chr_translocation)  # Remove helper column
    
    final_count <- nrow(igcaller_tsv)
    
    # Report filtering statistics
    message("  IgCaller filtering summary:")
    message("    Original: ", original_count, " SVs")
    message("    After base filter (Score>=10, ReadsNormal<=1, PoN<=1): ", base_filtered_count, " SVs")
    message("    After sex chr TRA filter (Score>=15, ReadsNormal=0, PoN=0): ", final_count, " SVs")
    message("    Removed: ", original_count - final_count, " SVs")
    
    result <- process_caller("igcaller", igcaller_tsv, "tsv", tum_norm_id, ref_seq_info, temp_files)
    if (result$has_data) {
      caller_jnc_list[["igcaller"]] <- result$jnc
      available_callers <- c(available_callers, "igcaller")
    }
    temp_files <- result$temp_files
  }
  
  # Check if we have at least one caller
  if (length(available_callers) == 0) {
    message("\n  WARNING: No SVs found across all callers for sample ", sample_set$sample[i], ". Skipping...")
    next
  }
  
  message("\nFound SVs in ", length(available_callers), " caller(s): ", paste(available_callers, collapse = ", "))
  
  # Dynamic consensus merge
  message("\nPerforming coordinate-based merge to create consensus junction set ...")
  consensus_jnc <- suppressWarnings(merge_callers_dynamic(caller_jnc_list, available_callers))
  
  # Convert gGnome junctions to BEDPE
  message("\nWriting output consensus gGnome junctions as BEDPE format ...")
  
  # Convert GRangeList of junctions to GRanges
  consensus_jnc_gr <- gUtils::grl.unlist(grl = consensus_jnc$grl)
  consensus_jnc_gr <- gr_refactor_seqs(input_gr = consensus_jnc_gr)
  
  # Convert GRanges to dataframe
  consensus_jnc_dt <- gUtils::gr2dt(consensus_jnc_gr)

  # Debug: Check the structure
  message("  ", nrow(consensus_jnc_dt) / 2, " BEDPE records")
  #message("  Debug: Column names: ", paste(colnames(consensus_jnc_dt), collapse = ", "))
  
  # Check if we have any merged junctions
  if (is.null(consensus_jnc_dt$merged.ix) || length(unique(consensus_jnc_dt$merged.ix)) == 0) {
    message("\n  WARNING: No merged junctions found. Check input data.")
    next
  }

  # Construct the final BEDPE output
  final_bedpe <- data.table::data.table()
  
  for (j in 1:length(unique(consensus_jnc_dt$merged.ix))) {
    
    # Get breakpoints of SV pair
    bedpe_pair <- consensus_jnc_dt %>%
      dplyr::filter(merged.ix == unique(consensus_jnc_dt$merged.ix)[j]) 

    # Debug check
    if (nrow(bedpe_pair) != 2) {
      message("  Warning: Expected 2 breakpoints for merged.ix ", unique(consensus_jnc_dt$merged.ix)[j], 
              " but found ", nrow(bedpe_pair))
      next
    }

    # Select relevant columns
    bedpe_pair <- bedpe_pair %>%
      dplyr::select(seqnames, start, end, strand,
                   dplyr::all_of(which(stringr::str_detect(string = colnames(consensus_jnc_dt), pattern = "name\\."))),
                   dplyr::all_of(which(stringr::str_detect(string = colnames(consensus_jnc_dt), pattern = "seen\\.by\\."))))
    
    # Find index of dynamic per caller columns
    name_col_idx <- which(stringr::str_detect(string = colnames(bedpe_pair), pattern = "name\\."))
    
    # Collect type of SV junction
    if (length(name_col_idx) > 0) {
      record_type <- stats::na.omit(bedpe_pair %>% 
                                   dplyr::select(dplyr::all_of(name_col_idx)) %>% 
                                   unique() %>% 
                                   t() %>% 
                                   as.vector())[1]
    } else {
      record_type <- "BND"  # Default if no type found
    }
    
    # Extract caller consensus string dynamically
    record_caller_combo_string <- extract_caller_string(bedpe_pair, available_callers)
    
    # Build the final BEDPE output
    bedpe_single_line_record <- data.table::data.table(
      "#chr1" = as.character(bedpe_pair$seqnames[1]),
      "start1" = bedpe_pair$start[1],
      "end1" = bedpe_pair$end[1],
      "chr2" = as.character(bedpe_pair$seqnames[2]),
      "start2" = bedpe_pair$start[2],
      "end2" = bedpe_pair$end[2],
      "type" = as.character(record_type),
      "score" = ".",
      "strand1" = as.character(bedpe_pair$strand[2]),
      "strand2" = as.character(bedpe_pair$strand[1]),
      "caller" = record_caller_combo_string
    )
    
    final_bedpe <- rbind(final_bedpe, bedpe_single_line_record)
  }

  # Check if we got any results
  if (nrow(final_bedpe) == 0) {
    message("\n  WARNING: No BEDPE records generated for sample ", sample_set$sample[i])
    next
  }
  
  # Write final BEDPE output file
  polished_bedpe_filename <- paste0(output_dir, tum_norm_id, ".hq.union.consensus.somatic.sv.bedpe")
  data.table::fwrite(x = final_bedpe,
                    file = polished_bedpe_filename,
                    sep = "\t",
                    col.names = TRUE,
                    quote = FALSE,
                    row.names = FALSE)
  
  message("\n", basename(polished_bedpe_filename), " ..... D O N E")
  
  # Clean up temp files
  unlink(temp_files)
}

