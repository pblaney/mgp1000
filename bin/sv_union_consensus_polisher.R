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

# Read in VCFs
message("Reading in SV calls.....")
manta_vcf <- VariantAnnotation::readVcf(file = manta_vcf_path, genome = "hg38")
delly_vcf <- VariantAnnotation::readVcf(file = delly_vcf_path, genome = "hg38")
svaba_vcf <- VariantAnnotation::readVcf(file = svaba_vcf_path, genome = "hg38")

# Adjust SvBAB VCF to include END info field
for(i in 1:nrow(VariantAnnotation::info(svaba_vcf))) {
  
  if(svaba_vcf@info$SPAN[i] != -1) {
    
    end_of_sv <- svaba_vcf@rowRanges[i]@ranges@start + svaba_vcf@info$SPAN[i]
    svaba_vcf@info$END[i] <- end_of_sv
    
  } else if(svaba_vcf@info$SPAN[i] == -1) {
    svaba_vcf@info$END[i] <- NA
  }
}

# Convert to breakpointRanges object, don't adjust for CIPOS uncertainty (i.e. keep nominalPosition)
manta_bp_ranges <- StructuralVariantAnnotation::breakpointRanges(manta_vcf, nominalPosition = T, inferMissingBreakends = TRUE)
svaba_bp_ranges <- StructuralVariantAnnotation::breakpointRanges(svaba_vcf, nominalPosition = T, inferMissingBreakends = TRUE)
delly_bp_ranges <- StructuralVariantAnnotation::breakpointRanges(delly_vcf, nominalPosition = T, inferMissingBreakends = TRUE)

# Read in IgCaller oncogeneic rearrangements TSV, filter then convert to BEDPE
igcaller <- readr::read_delim(file = igcaller_tsv_path,
                              delim = "\t",
                              col_names = T,
                              show_col_types = F)

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
                                  "strand2" = StrandB)

igcaller_bedpe$type <- dplyr::case_when(igcaller_bedpe$type == "Deletion" ~ "DEL",
                                        igcaller_bedpe$type == "Duplication" ~ "DUP",
                                        igcaller_bedpe$type == "Insertion" ~ "INS",
                                        igcaller_bedpe$type == "Inversion" ~ "INV",
                                        igcaller_bedpe$type == "Translocation" ~ "TRA")

# Convert the breakpointRanges object to BEDPE format
message("Converting SV calls to common BEDPE format.....")
manta_bedpe <- vcfToBedpe(vcf = manta_bp_ranges)
svaba_bedpe <- vcfToBedpe(vcf = svaba_bp_ranges)
delly_bedpe <- vcfToBedpe(vcf = delly_bp_ranges)

# Write out temporary BEDPE files
temp_manta_bedpe_filename <- tempfile(pattern = tum_norm_id, fileext = ".manta.somatic.sv.temp.bedpe")
write.table(x = manta_bedpe, file = temp_manta_bedpe_filename, row.names = F, col.names = T, sep = "\t", quote = F)

temp_delly_bedpe_filename <- tempfile(pattern = tum_norm_id, fileext = ".delly.somatic.sv.temp.bedpe")
write.table(x = delly_bedpe, file = temp_delly_bedpe_filename, row.names = F, col.names = T, sep = "\t", quote = F)

temp_svaba_bedpe_filename <- tempfile(pattern = tum_norm_id, fileext = ".svaba.somatic.sv.temp.bedpe")
write.table(x = svaba_bedpe, file = temp_svaba_bedpe_filename, row.names = F, col.names = T, sep = "\t", quote = F)

temp_igcaller_bedpe_filename <- tempfile(pattern = tum_norm_id, fileext = ".igcaller.somatic.sv.temp.bedpe")
write.table(x = igcaller_bedpe, file = temp_igcaller_bedpe_filename, row.names = F, col.names = T, sep = "\t", quote = F)

# Read in BEDPE files into gGnome junction objects
message("Converting BEDPE to consensus gGnome junctions.....")
ref_seq_info <- c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
                  159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
                  114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
                  58617616, 64444167, 46709983, 50818468, 156040895, 57227415)
names(ref_seq_info) <- paste0("chr", c(seq(1,22,1), "X", "Y"))

manta_jnc <- gGnome::jJ(rafile = temp_manta_bedpe_filename,
                		chr.convert = F,
                		hg = "hg38",
                		keep.features = T,
                		seqlengths = ref_seq_info)

delly_jnc <- gGnome::jJ(rafile = temp_delly_bedpe_filename,
                		chr.convert = F,
                		hg = "hg38",
                		keep.features = T,
                		seqlengths = ref_seq_info)

svaba_jnc <- gGnome::jJ(rafile = temp_svaba_bedpe_filename,
                		chr.convert = F,
                		hg = "hg38",
                		keep.features = T,
                		seqlengths = ref_seq_info)

igcaller_jnc <- gGnome::jJ(rafile = temp_igcaller_bedpe_filename,
                  		   chr.convert = F,
                  		   hg = "hg38",
                  		   keep.features = T,
                  		   seqlengths = ref_seq_info)

# Comprehensive coordinate-based merge to create consensus junction set
consensus_jnc <- gGnome::merge(manta = manta_jnc[,"name"],
							                 delly = delly_jnc[,"name"],
                               svaba = svaba_jnc[,"name"],
                               igcaller = igcaller_jnc[,"name"],
							                 pad = 1000)

# Generate intersection UpSetR plot
plot_data <- as.data.frame(sign(as.matrix(consensus_jnc$dt[,.(seen.by.manta, seen.by.delly, seen.by.svaba, seen.by.igcaller)])))
pdf(file = paste0(tum_norm_id, ".hq.union.consensus.somatic.sv.intersection.plot.pdf"),
    width = 7,
    height = 7)
UpSetR::upset(plot_data,
              mainbar.y.label = "Breakpoint Pair Intersection",
              sets.x.label = "SVs Per Caller",
              point.size = 5.25,
              line.size = 1.8,
              text.scale = c(1.75, 1.5, 1.5, 1.25, 1.25, 1.5),
              scale.intersections = "identity",
              scale.sets = "identity",
              main.bar.color = paletteer::paletteer_dynamic("cartography::wine.pal",
                                                            n = length(table(paste0(plot_data$seen.by.manta,
                                                                                    plot_data$seen.by.delly,
                                                                                    plot_data$seen.by.svaba,
                                                                                    plot_data$seen.by.igcaller)))),
              sets.bar.color = paletteer::paletteer_d("ggsci::default_aaas")[3:6]) # BigY, Yaxistix, LittleX, Xaxistix, callerId, colcountannots
dev.off()
message("Writing output consensus gGnome junctions as BEDPE format.....")
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
                  name.manta, name.delly, name.svaba, name.igcaller,
                  seen.by.manta, seen.by.delly, seen.by.svaba, seen.by.igcaller)
  
  record_type <- na.omit(unique(c(bedpe_pair$name.manta, bedpe_pair$name.delly, bedpe_pair$name.svaba, bedpe_pair$name.igcaller)))[1]
  
  # Collect evaluation of seen.by columns to generate a string of caller consensus
  record_caller <- data.table::data.table(manta = bedpe_pair$seen.by.manta,
                                          delly = bedpe_pair$seen.by.delly,
                                          svaba = bedpe_pair$seen.by.svaba,
                                          igcaller = bedpe_pair$seen.by.igcaller)
  record_caller <- colSums(record_caller)
  
  record_caller_string <- c(NA, NA, NA, NA)
  record_caller_string[1] <- case_when(record_caller["manta"] > 0 ~ "manta", .default = NA)
  record_caller_string[2] <- case_when(record_caller["delly"] > 0 ~ "delly", .default = NA)
  record_caller_string[3] <- case_when(record_caller["svaba"] > 0 ~ "svaba", .default = NA)
  record_caller_string[4] <- case_when(record_caller["igcaller"] > 0 ~ "igcaller", .default = NA)
  
  if(length(which(!is.na(record_caller_string))) > 1) {
    record_caller_combo_string <- stringr::str_flatten(string = record_caller_string[which(!is.na(record_caller_string))], collapse = ",")
  } else {
    record_caller_combo_string <- record_caller_string[which(!is.na(record_caller_string))]
  }
  
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

unlink(temp_manta_bedpe_filename)
unlink(temp_delly_bedpe_filename)
unlink(temp_svaba_bedpe_filename)
unlink(temp_igcaller_bedpe_filename)

