#!/usr/bin/env Rscript
# This script was original writen by Alex Di Genova at https://github.com/IARCbioinfo/facets-nf/blob/master/bin/facets.cval.r
# It has been adapted here

library(facets)
library(data.table)
options(datatable.fread.input.cmd.message = FALSE)
source(system.file("extRfns", "readSnpMatrixDT.R", package = "facets"))

datafile = commandArgs(TRUE)[1] # "MESO_054_T.chr1.pileup.gz"
genome = commandArgs(TRUE)[2] # "hg19"
cur_params = as.numeric(c(commandArgs(TRUE)[3:6])) # c(1000,25,300,150)
ndepth_param = as.numeric(commandArgs(TRUE)[7]) # 20

# Will provide multiple cvalues and output PDF by default
m_cval = TRUE
plot_pdf = TRUE

# Check if we have to compute multiples cvalues
#if (!is.na(commandArgs(TRUE)[8]) && commandArgs(TRUE)[8]=="MCVAL") {
#	m_cval = TRUE
#} else {#
#	m_cval = FALSE
#}

#if (!is.na(commandArgs(TRUE)[9]) && commandArgs(TRUE)[9]=="PDF") {
#	plot_pdf = TRUE
#} else {
#	plot_pdf = FALSE
#}

#we set a random set for reproducible analysis PreprocSamples
set.seed(2022)

sample_name = gsub("csv.gz", "", datafile)
sample_name_no_dot = gsub(".csv.gz", "", datafile)
rcmat = readSnpMatrixDT(datafile)

plot_facets = function (oo_facets, fit_facets, text_title, plot_name, pref = "cval300", pdf = F) {
  if(pdf) {
  	pdf(paste(plot_name,pref, "_CNV.pdf", sep = ""), width = 12, height = 11)
  } else {
  	png(paste(plot_name,pref, "_CNV.png", sep = ""), width = 30, height = 27.5, units = "cm", res = 300)
  }
  plotSample(x = oo_facets, emfit = fit_facets, sname = text_title)
  dev.off()
  pdf(paste(plot_name, pref, "_CNV_spider.pdf", sep = ""), width = 6, height = 6)
  logRlogORspider(oo_facets$out, oo_facets$dipLogR)
  dev.off()
}


# Default parameters for typical +30x coverage tumor: snp_nbhd = 1000, cval_preproc = 35, cval_proc1 = 300, cval_proc2 = 150, min_read_count = 20
# Note, snp_nbhd should range between 500-5000, where a 500 window size would lead to the most number of SNPs used
# cval_preproc :
xx = preProcSample(rcmat, gbuild = genome, ndepth = ndepth_param, snp.nbhd = cur_params[1], cval = cur_params[2])
fit1=procSample(xx, cval = cur_params[3], min.nhet = 15, dipLogR=NULL) # fit fine

# You can specify cval that is large enough to avoid hyper-segmentation
oo_fine = procSample(xx, cval = cur_params[4], min.nhet = 15, dipLogR = fit1$dipLogR) # large fit to avoid short segments
fit_fine = emcncf(oo_fine)

# Plot the results
text_title=paste(sample_name_no_dot, ": Purity=", round(fit_fine$purity, 3) * 100, "%; Ploidy=", round(fit_fine$ploidy, 2), sep = "")
pref=paste("def_cval", cur_params[4], sep = "")
plot_facets(oo_fine, fit_fine, text_title, sample_name, pref, plot_pdf)

# Write the output table
cat("", "purity", "ploidy", "dipLogR", "loglik", "\n", file = paste(sample_name, pref, "_stats.txt", sep = ""), sep= "\t")
cat(sample_name_no_dot, fit_fine$purity, fit_fine$ploidy, fit_fine$dipLogR, fit_fine$loglik, file = paste(sample_name, pref, "_stats.txt", sep = ""), sep = "\t", append = T)
fit_fine$cncf['cnlr.median-dipLogR'] = fit_fine$cncf$cnlr.median - fit_fine$dipLogR
write.table(fit_fine$cncf, file = paste(sample_name, pref, "_CNV.txt", sep = ""), quote = F, sep = "\t", row.names = F)

# Compute multiples cvalues 500, 1000, and 1500
if(m_cval) {
	
  # You can specify cval that is large enough to avoid hyper-segmentation
  oo_fine = procSample(xx, cval = 500, min.nhet = 15, dipLogR = fit1$dipLogR) # large fit to avoid short segments
  fit_fine = emcncf(oo_fine)
  
  # Plot the results
  text_title = paste(sample_name_no_dot, ": Purity=", round(fit_fine$purity, 3) * 100, "%; Ploidy=", round(fit_fine$ploidy, 2), sep = "")
  plot_facets(oo_fine, fit_fine, text_title, sample_name, "cval500", plot_pdf)

  # Write the ouput table
  cat("", "purity", "ploidy", "dipLogR", "loglik", "\n", file = paste(sample_name, "cval500", "_stats.txt", sep = ""), sep= "\t")
  cat(sample_name_no_dot, fit_fine$purity, fit_fine$ploidy, fit_fine$dipLogR, fit_fine$loglik, file = paste(sample_name, "cval500", "_stats.txt", sep = ""), sep = "\t", append = T)
  fit_fine$cncf['cnlr.median-dipLogR'] = fit_fine$cncf$cnlr.median - fit_fine$dipLogR
  write.table(fit_fine$cncf, file = paste(sample_name, "cval500", "_CNV.txt", sep = ""), quote = F, sep = "\t", row.names = F)

  # You can specify cval that is large enough to avoid hyper-segmentation
  oo_fine = procSample(xx, cval = 1000, min.nhet = 15, dipLogR = fit1$dipLogR) # large fit to avoid short segments
  fit_fine = emcncf(oo_fine)

  # Plot the results
  text_title = paste(sample_name_no_dot, ": Purity=", round(fit_fine$purity, 3) * 100, "%; Ploidy=", round(fit_fine$ploidy, 2), sep = "")
  plot_facets(oo_fine, fit_fine, text_title, sample_name, "cval1000", plot_pdf)

  # Write the ouput table
  cat("", "purity", "ploidy", "dipLogR", "loglik", "\n", file = paste(sample_name, "cval1000", "_stats.txt", sep = ""), sep = "\t")
  cat(sample_name_no_dot, fit_fine$purity, fit_fine$ploidy, fit_fine$dipLogR, fit_fine$loglik, file = paste(sample_name, "cval1000", "_stats.txt", sep = ""), sep = "\t", append = T)
  fit_fine$cncf['cnlr.median-dipLogR'] = fit_fine$cncf$cnlr.median - fit_fine$dipLogR
  write.table(fit_fine$cncf, file = paste(sample_name, "cval1000", "_CNV.txt", sep = ""), quote = F, sep = "\t", row.names = F)

  # You can specify cval that is large enough to avoid hyper-segmentation
  oo_fine=procSample(xx, cval = 1500, min.nhet = 15, dipLogR = fit1$dipLogR) # large fit to avoid short segments
  fit_fine = emcncf(oo_fine)

  # Plot the results
  text_title = paste(sample_name_no_dot, ": Purity=", round(fit_fine$purity, 3) * 100, "%; Ploidy=", round(fit_fine$ploidy, 2), sep = "")
  plot_facets(oo_fine, fit_fine, text_title, sample_name, "cval1500", plot_pdf)

  #Write the ouput table
  cat("", "purity", "ploidy", "dipLogR", "loglik", "\n", file = paste(sample_name, "cval1500", "_stats.txt", sep = ""), sep = "\t")
  cat(sample_name_no_dot, fit_fine$purity, fit_fine$ploidy, fit_fine$dipLogR, fit_fine$loglik, file = paste(sample_name, "cval1500", "_stats.txt", sep = ""), sep = "\t", append = T)
  fit_fine$cncf['cnlr.median-dipLogR'] = fit_fine$cncf$cnlr.median - fit_fine$dipLogR
  write.table(fit_fine$cncf, file = paste(sample_name, "cval1500", "_CNV.txt", sep = ""), quote = F, sep = "\t", row.names = F)

}

# Save the facets info
writeLines(capture.output(sessionInfo()), paste(sample_name, "R_sessionInfo.txt", sep = ""))
