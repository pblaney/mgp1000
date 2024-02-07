#!/usr/bin/env Rscript

# This script accepts set of SV VCF files as input and generates
# a union consensus SV file

#########################
#####   Libraries   #####

suppressPackageStartupMessages(library(StructuralVariantAnnotation))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(gGnome))
suppressPackageStartupMessages(library(UpSetR))
suppressPackageStartupMessages(library(devgru))
suppressPackageStartupMessages(library(paletteer))
suppressPackageStartupMessages(library(ggsci))

Sys.setenv(DEFAULT_BSGENOME = 'BSgenome.Hsapiens.UCSC.hg38::Hsapiens')
options(scipen = 999)

#########################
#####   Functions   #####

# TODO: Convert this to function that wraps all parts together to go from VCF to BEDPE
# Original author Jennifer Shelton, https://bitbucket.nygenome.org/projects/WDL/repos/somatic_dna_tools/browse/vcf_to_bedpe.r
# commit ef47c0317e8  14 Jan 2022
## Convert breakpointRanges to BEDPE
vcfToBedpe = function(vcf) {
  
  sqn = as.character(seqnames(vcf))
  strand = as.character(strand(vcf))
  res = c()
  processed = c()
  
  
  for (i in 1:length(vcf)) {
    bnd = names(vcf)[i]
    partner = vcf$partner[i]
    partner.idx = which(names(vcf) == partner)
    
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
    res.i = c(sqn[i], start(vcf)[i], end(vcf)[i],                                  ## chr1, start1, end1
              sqn[partner.idx], start(vcf)[partner.idx], end(vcf)[partner.idx],    ## chr2, start2, end 2
              'BND', '.', strand[i], strand[partner.idx])                 ## type, score, strand1, strand2, support
    
    ## Add to result, keep track of processed breakends
    res = rbind(res, res.i)
    processed = c(processed, bnd, partner)
  }
  
  ## Add colnames and fill in simple event classifications
  colnames(res) = c('chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'type', 'score', 'strand1', 'strand2')
  res = as.data.frame(res, stringsAsFactors=F)
  
  res$type[res$strand1 == '+' & res$strand2 == '-'] = 'DEL'
  res$type[res$strand1 == '-' & res$strand2 == '+'] = 'DUP'
  res$type[res$strand1 == '-' & res$strand2 == '-'] = 'INV'
  res$type[res$strand1 == '+' & res$strand2 == '+'] = 'INV'
  res$type[res$chr1 != res$chr2] = 'TRA'
  
  ## Sort by chromosome 
  res = res[order(factor(res$chr1, levels=levels(seqnames(vcf))), res$start1, res$end1, decreasing=F), ]
  
  ## Simplify coordinates
  res$end1 = as.numeric(res$start1) + 1
  res$end2 = as.numeric(res$start2) + 1
  
  
  colnames(res)[1] = paste0('#', colnames(res)[1])
  
  return(res)
  
}

#########################
#####   Execution   #####

# Accept command line arguments as input
input_args <- commandArgs(trailingOnly = T)
  
tum_norm_id <- input_args[1]

manta_vcf_path <- input_args[2]

svaba_vcf_path <- input_args[3]

delly_vcf_path <- input_args[4]

igcaller_tsv_path <- input_args[5]

# Set seqinfo named vector of hg38 chrom names and lengths
ref_seq_info <- c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
                  159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
                  114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
                  58617616, 64444167, 46709983, 50818468, 156040895, 57227415)
names(ref_seq_info) <- paste0("chr", c(seq(1,22,1), "X", "Y"))

# Read in SV calls
message("Reading in SV calls ...")
manta_vcf <- VariantAnnotation::readVcf(file = manta_vcf_path,
                                        genome = "hg38")
svaba_vcf <- VariantAnnotation::readVcf(file = svaba_vcf_path,
                                        genome = "hg38")
delly_vcf <- VariantAnnotation::readVcf(file = delly_vcf_path,
                                        genome = "hg38")
igcaller_tsv <- readr::read_delim(file = igcaller_tsv_path,
                                  delim = "\t",
                                  col_names = T,
                                  show_col_types = F)

# Filter IgCaller calls
igcaller_tsv <- igcaller_tsv %>%
                  dplyr::filter(Score >= 5 & `Reads in normal` <= 4 & `Count in PoN` <= 5)

# For each input file, first check for empty files where there are SVs surviving caller specific filtering (seen in Manta and IgCaller),
# If not empty, convert each caller SV set to breakpointRanges object, don't adjust for CIPOS uncertainty (i.e. keep nominalPosition)
# then convert the breakpointRanges object to BEDPE format and finally convert these to gGnome junction objects
processed_input <- ""
temp_files <- c()
message("Processing Manta input ...")
message("  Found ", nrow(VariantAnnotation::info(manta_vcf)), " SVs in VCF ...")
if(nrow(VariantAnnotation::info(manta_vcf)) > 0) {
  
  # VCF ---> breakpointRanges
  manta_bp_ranges <- StructuralVariantAnnotation::breakpointRanges(manta_vcf,
                                                                   nominalPosition = T,
                                                                   inferMissingBreakends = TRUE)
  # breakpointRanges ---> BEDPE
  manta_bedpe <- vcfToBedpe(vcf = manta_bp_ranges)
  temp_manta_bedpe_filename <- tempfile(pattern = tum_norm_id,
                                        fileext = ".manta.somatic.sv.temp.bedpe")
  # add temp file to tracker and write out temp BEDPE
  temp_files <- append(temp_files, temp_manta_bedpe_filename)
  write.table(x = manta_bedpe, file = temp_manta_bedpe_filename,
              row.names = F,
              col.names = T,
              sep = "\t",
              quote = F)
  # BEDPE ---> gGnome Juncitons
  manta_jnc <- gGnome::jJ(rafile = temp_manta_bedpe_filename,
                          chr.convert = F,
                          hg = "hg38",
                          keep.features = T,
                          seqlengths = ref_seq_info)
  
  # Track if processed
  processed_input <- "manta,"
  
} else {
  # If no records, set to NA
  manta_jnc <- NA
}

message("Processing SvABA input ...")
message("  Found ", nrow(VariantAnnotation::info(svaba_vcf)), " SVs in VCF ...")
if(nrow(VariantAnnotation::info(svaba_vcf)) > 0) {
  
  # Adjust SvBAB VCF to include END info field
  for(i in 1:nrow(VariantAnnotation::info(svaba_vcf))) {
    
    if(svaba_vcf@info$SPAN[i] != -1) {
      end_of_sv <- svaba_vcf@rowRanges[i]@ranges@start + svaba_vcf@info$SPAN[i]
      svaba_vcf@info$END[i] <- end_of_sv
      
    } else if(svaba_vcf@info$SPAN[i] == -1) {
      svaba_vcf@info$END[i] <- NA
    }
  }
  
  # VCF ---> breakpointRanges
  svaba_bp_ranges <- StructuralVariantAnnotation::breakpointRanges(svaba_vcf,
                                                                   nominalPosition = T,
                                                                   inferMissingBreakends = TRUE)
  # breakpointRanges ---> BEDPE
  svaba_bedpe <- vcfToBedpe(vcf = svaba_bp_ranges)
  temp_svaba_bedpe_filename <- tempfile(pattern = tum_norm_id,
                                        fileext = ".svaba.somatic.sv.temp.bedpe")
  # add temp file to tracker and write out temp BEDPE
  temp_files <- append(temp_files, temp_svaba_bedpe_filename)
  write.table(x = svaba_bedpe,
              file = temp_svaba_bedpe_filename,
              row.names = F,
              col.names = T,
              sep = "\t",
              quote = F)
  # BEDPE ---> gGnome Juncitons
  svaba_jnc <- gGnome::jJ(rafile = temp_svaba_bedpe_filename,
                          chr.convert = F,
                          hg = "hg38",
                          keep.features = T,
                          seqlengths = ref_seq_info)
  
  # Track if processed
  processed_input <- stringr::str_c(processed_input, "svaba,")
  
} else {
  # If no records, set to NA
  svaba_jnc <- NA
}

message("Processing DELLY2 input ...")
message("  Found ", nrow(VariantAnnotation::info(delly_vcf)), " SVs in VCF ...")
if(nrow(VariantAnnotation::info(delly_vcf)) > 0) {
  
  # VCF ---> breakpointRanges
  delly_bp_ranges <- StructuralVariantAnnotation::breakpointRanges(delly_vcf,
                                                                   nominalPosition = T,
                                                                   inferMissingBreakends = TRUE)
  # breakpointRanges ---> BEDPE
  delly_bedpe <- vcfToBedpe(vcf = delly_bp_ranges)
  temp_delly_bedpe_filename <- tempfile(pattern = tum_norm_id,
                                        fileext = ".delly.somatic.sv.temp.bedpe")
  # add temp file to tracker and write out temp BEDPE
  temp_files <- append(temp_files, temp_delly_bedpe_filename)
  write.table(x = delly_bedpe,
              file = temp_delly_bedpe_filename,
              row.names = F,
              col.names = T,
              sep = "\t",
              quote = F)
  # BEDPE ---> gGnome Juncitons
  delly_jnc <- gGnome::jJ(rafile = temp_delly_bedpe_filename,
                          chr.convert = F,
                          hg = "hg38",
                          keep.features = T,
                          seqlengths = ref_seq_info)
  
  # Track if processed
  processed_input <- stringr::str_c(processed_input, "delly,")
  
} else {
  # If no records, set to NA
  delly_jnc <- NA
}

message("Processing IgCaller input ...")
message("  Found ", nrow(igcaller_tsv), " SVs in TSV ...")
if(nrow(igcaller_tsv) > 0) {
  
  # TSV ---> BEDPE
  igcaller_bedpe <- igcaller_tsv %>%
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

  igcaller_bedpe$type <- dplyr::case_when(igcaller_bedpe$type == "Deletion" ~ "DEL",
                                          igcaller_bedpe$type == "Duplication" ~ "DUP",
                                          igcaller_bedpe$type == "Gain" ~ "DUP",
                                          igcaller_bedpe$type == "Insertion" ~ "INS",
                                          igcaller_bedpe$type == "Inversion" ~ "INV",
                                          igcaller_bedpe$type == "Translocation" ~ "TRA")
  
  temp_igcaller_bedpe_filename <- tempfile(pattern = tum_norm_id,
                                           fileext = ".igcaller.somatic.sv.temp.bedpe")
  # add temp file to tracker and write out temp BEDPE
  temp_files <- append(temp_files, temp_igcaller_bedpe_filename)
  write.table(x = igcaller_bedpe,
              file = temp_igcaller_bedpe_filename,
              row.names = F,
              col.names = T,
              sep = "\t",
              quote = F)
  # BEDPE ---> gGnome Juncitons
  igcaller_jnc <- gGnome::jJ(rafile = temp_igcaller_bedpe_filename,
                         chr.convert = F,
                         hg = "hg38",
                         keep.features = T,
                         seqlengths = ref_seq_info)
  
  # Track if processed
  processed_input <- stringr::str_c(processed_input, "igcaller")
  
} else {
  # If no records, set to NA
  igcaller_jnc <- NA
}

# Comprehensive coordinate-based merge to create consensus junction set
if(processed_input == "manta,svaba,delly,igcaller") {
  # Comprehensive coordinate-based merge to create consensus junction set
  consensus_jnc <- gGnome::merge(manta = manta_jnc[,"name"],
                                 svaba = svaba_jnc[,"name"],
                                 delly = delly_jnc[,"name"],
                                 igcaller = igcaller_jnc[,"name"],
                                 pad = 1000)
  # Set plot data for UpSetR intersection plot
  plot_data <- as.data.frame(sign(as.matrix(consensus_jnc$dt[,.(seen.by.manta, seen.by.svaba, seen.by.delly, seen.by.igcaller)])))
  plot_bar_cols <- length(table(paste0(plot_data$seen.by.manta,
                                       plot_data$seen.by.svaba,
                                       plot_data$seen.by.delly,
                                       plot_data$seen.by.igcaller)))
  
  # Edge Case: No surviving records in Manta (even after base filtering, i.e. bland genomes)
} else if(processed_input == "svaba,delly,igcaller") {
  # Comprehensive coordinate-based merge to create consensus junction set
  consensus_jnc <- gGnome::merge(svaba = svaba_jnc[,"name"],
                                 delly = delly_jnc[,"name"],
                                 igcaller = igcaller_jnc[,"name"],
                                 pad = 1000)
  # Set plot data for UpSetR intersection plot
  plot_data <- as.data.frame(sign(as.matrix(consensus_jnc$dt[,.(seen.by.svaba, seen.by.delly, seen.by.igcaller)])))
  plot_bar_cols <- length(table(paste0(plot_data$seen.by.svaba,
                                       plot_data$seen.by.delly,
                                       plot_data$seen.by.igcaller)))
  
  # Edge Case: No surviving records in IgCaller (even after base filtering, i.e. bland genomes)
} else if(processed_input == "manta,svaba,delly,") {
  # Comprehensive coordinate-based merge to create consensus junction set
  consensus_jnc <- gGnome::merge(manta = manta_jnc[,"name"],
                                 svaba = svaba_jnc[,"name"],
                                 delly = delly_jnc[,"name"],
                                 pad = 1000)
  # Set plot data for UpSetR intersection plot
  plot_data <- as.data.frame(sign(as.matrix(consensus_jnc$dt[,.(seen.by.manta, seen.by.svaba, seen.by.delly)])))
  plot_bar_cols <- length(table(paste0(plot_data$seen.by.manta,
                                       plot_data$seen.by.svaba,
                                       plot_data$seen.by.delly)))
  
  # Edge Case: No surviving records in Manta/IgCaller
} else if(processed_input == "svaba,delly,") {
  # Comprehensive coordinate-based merge to create consensus junction set
  consensus_jnc <- gGnome::merge(svaba = svaba_jnc[,"name"],
                                 delly = delly_jnc[,"name"],
                                 pad = 1000)
  # Set plot data for UpSetR intersection plot
  plot_data <- as.data.frame(sign(as.matrix(consensus_jnc$dt[,.(seen.by.svaba, seen.by.delly)])))
  plot_bar_cols <- length(table(paste0(plot_data$seen.by.svaba,
                                       plot_data$seen.by.delly)))
}

message("Generating UpSet plot of SV breakpoint pair intersection ...")
pdf(file = paste0(tum_norm_id, ".hq.union.consensus.somatic.sv.intersection.plot.pdf"),
    width = 7,
    height = 7,
    onefile = FALSE)
UpSetR::upset(plot_data,
              mainbar.y.label = "Breakpoint Pair Intersection",
              sets.x.label = "SVs Per Caller",
              point.size = 5.25,
              line.size = 1.8,
              text.scale = c(1.75, 1.5, 1.5, 1.25, 1.25, 1.5),
              scale.intersections = "identity",
              scale.sets = "identity",
              sets.bar.color = paletteer::paletteer_d("ggsci::default_aaas")[3:(2 + length(which(list(manta_jnc, svaba_jnc, delly_jnc, igcaller_jnc) != "NA")))],
              main.bar.color = paletteer::paletteer_dynamic("cartography::wine.pal",
                                                            n = plot_bar_cols))
dev.off()

message("Writing output consensus gGnome junctions as BEDPE format ...")
# TODO: Go from SV gGnome junctions to BEDPE
# Convert GRangeList of junctions to GRanges
consensus_jnc_gr <- gUtils::grl.unlist(grl = consensus_jnc$grl)
consensus_jnc_gr <- gr_refactor_seqs(input_gr = consensus_jnc_gr)

# Convert GRanges to dataframe for easy manipulation
consensus_jnc_dt <- gUtils::gr2dt(consensus_jnc_gr)

# Construct the final BEDPE output
final_bedpe <- tibble::tibble()
for (i in 1:length(unique(consensus_jnc_dt$merged.ix))) {
  
  # Get breakpoints of SV pair
  bedpe_pair <- consensus_jnc_dt %>%
    dplyr::filter(merged.ix == unique(consensus_jnc_dt$merged.ix)[i]) %>%
    dplyr::select(seqnames, start, end, strand,
                  dplyr::all_of(which(stringr::str_detect(string = colnames(consensus_jnc_dt), pattern = "name\\."))),
                  dplyr::all_of(which(stringr::str_detect(string = colnames(consensus_jnc_dt), pattern = "seen\\.by\\."))))
  
  # Find index of dynamic per caller columns
  seen_by_col_idx <- which(stringr::str_detect(string = colnames(bedpe_pair), pattern = "seen\\.by\\."))
  name_col_idx <- which(stringr::str_detect(string = colnames(bedpe_pair), pattern = "name\\."))
  
  # Collect type of SV junction and evaluation of seen.by columns to generate a string of caller consensus
  record_type <- stats::na.omit(bedpe_pair %>% dplyr::select(dplyr::all_of(name_col_idx)) %>% unique() %>% t() %>% as.vector())[1]
  
  # Comprehensive coordinate-based merge to create consensus junction set
  if(processed_input == "manta,svaba,delly,igcaller") {
    record_caller <- data.table::data.table(manta = bedpe_pair$seen.by.manta,
                                            svaba = bedpe_pair$seen.by.svaba,
                                            delly = bedpe_pair$seen.by.delly,
                                            igcaller = bedpe_pair$seen.by.igcaller)
    record_caller <- colSums(record_caller)
    
    record_caller_string <- c(NA, NA, NA, NA)
    record_caller_string[1] <- dplyr::case_when(record_caller["manta"] > 0 ~ "manta", .default = NA)
    record_caller_string[2] <- dplyr::case_when(record_caller["svaba"] > 0 ~ "svaba", .default = NA)
    record_caller_string[3] <- dplyr::case_when(record_caller["delly"] > 0 ~ "delly", .default = NA)
    record_caller_string[4] <- dplyr::case_when(record_caller["igcaller"] > 0 ~ "igcaller", .default = NA)
    
  } else if(processed_input == "svaba,delly,igcaller") {
    record_caller <- data.table::data.table(svaba = bedpe_pair$seen.by.svaba,
                                            delly = bedpe_pair$seen.by.delly,
                                            igcaller = bedpe_pair$seen.by.igcaller)
    record_caller <- colSums(record_caller)
    
    record_caller_string <- c(NA, NA, NA)
    record_caller_string[1] <- dplyr::case_when(record_caller["svaba"] > 0 ~ "svaba", .default = NA)
    record_caller_string[2] <- dplyr::case_when(record_caller["delly"] > 0 ~ "delly", .default = NA)
    record_caller_string[3] <- dplyr::case_when(record_caller["igcaller"] > 0 ~ "igcaller", .default = NA)
    
  } else if(processed_input == "manta,svaba,delly,") {
    record_caller <- data.table::data.table(manta = bedpe_pair$seen.by.manta,
                                            svaba = bedpe_pair$seen.by.svaba,
                                            delly = bedpe_pair$seen.by.delly)
    record_caller <- colSums(record_caller)
    
    record_caller_string <- c(NA, NA, NA)
    record_caller_string[1] <- dplyr::case_when(record_caller["manta"] > 0 ~ "manta", .default = NA)
    record_caller_string[2] <- dplyr::case_when(record_caller["svaba"] > 0 ~ "svaba", .default = NA)
    record_caller_string[3] <- dplyr::case_when(record_caller["delly"] > 0 ~ "delly", .default = NA)
    
  } else if(processed_input == "svaba,delly,") {
    
    record_caller <- data.table::data.table(svaba = bedpe_pair$seen.by.svaba,
                                            delly = bedpe_pair$seen.by.delly)
    record_caller <- colSums(record_caller)
    
    record_caller_string <- c(NA, NA, NA)
    record_caller_string[1] <- dplyr::case_when(record_caller["manta"] > 0 ~ "manta", .default = NA)
    record_caller_string[2] <- dplyr::case_when(record_caller["svaba"] > 0 ~ "svaba", .default = NA)
    record_caller_string[3] <- dplyr::case_when(record_caller["delly"] > 0 ~ "delly", .default = NA)
  }
  
  # Build the caller combo string for the final BEDPE output
  if(length(which(!is.na(record_caller_string))) > 1) {
    record_caller_combo_string <- stringr::str_flatten(string = record_caller_string[which(!is.na(record_caller_string))], collapse = ",")
  } else {
    record_caller_combo_string <- record_caller_string[which(!is.na(record_caller_string))]
  }
  
  # Build the final BEDPE output
  bedpe_single_line_record <- data.table::data.table("#chr1" = bedpe_pair$seqnames[1],
                                                     "start1" = bedpe_pair$start[1],
                                                     "end1" = bedpe_pair$end[1],
                                                     "chr2" = bedpe_pair$seqnames[2],
                                                     "start2" = bedpe_pair$start[2],
                                                     "end2" = bedpe_pair$end[2],
                                                     "type" = record_type,
                                                     "score" = ".",
                                                     "strand1" = bedpe_pair$strand[1],
                                                     "strand2" = bedpe_pair$strand[2],
                                                     "caller" = record_caller_combo_string)
  
  final_bedpe <- rbind(final_bedpe,
                       bedpe_single_line_record)
}

# Write final BEDPE output file
polished_bedpe_filename <- paste0(tum_norm_id, ".hq.union.consensus.somatic.sv.bedpe")

data.table::fwrite(x = final_bedpe,
                   file = polished_bedpe_filename,
                   sep = "\t",
                   col.names = TRUE,
                   quote = FALSE,
                   row.names = FALSE)

message(paste0(polished_bedpe_filename, "..... D O N E"))

unlink(temp_files)
