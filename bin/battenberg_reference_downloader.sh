#!/usr/bin/env bash
# Download all the necessary files for use in Battenberg

# Function
battenbergRefDownload() {

	# GC Correction
	wget -q -O GC_correction_hg38_chr.zip "https://www.dropbox.com/sh/bize1n830t0mgzb/AADQD4DTJOF75YmhBDDoQ9nla/GC_correction_hg38?dl=0&lst=" && \
  	unzip GC_correction_hg38_chr.zip
	mv 1000G_GC_chr*.txt.gz GC_correction_hg38/ && \
	rm GC_correction_hg38_chr.zip

	# RT Correction
  	wget -q -O RT_correction_hg38.zip "https://www.dropbox.com/sh/bize1n830t0mgzb/AABZ2uM13YMYB_q6X1pP1McJa/RT_correction_hg38?dl=0&lst=" && \
  	unzip RT_correction_hg38.zip
  	mv 1000G_RT_chr*.txt.gz RT_correction_hg38/ && \
  	rm RT_correction_hg38.zip

  	# Shapeit2
  	wget -q -O shapeit2_chr.zip "https://ora.ox.ac.uk/objects/uuid:08e24957-7e76-438a-bd38-66c48008cf52/download_file?file_format=&safe_filename=shapeit2_chr.zip&type_of_work=Dataset" && \
  	unzip shapeit2_chr.zip
  	rm shapeit2_chr.zip

  	# Impute
  	wget -q -O imputation_chr.zip "https://ora.ox.ac.uk/objects/uuid:08e24957-7e76-438a-bd38-66c48008cf52/download_file?file_format=&safe_filename=imputation_chr.zip&type_of_work=Dataset" && \
  	unzip imputation_chr.zip
  	rm imputation_chr.zip

  	# 1000G
  	wget -q -O 1000G_loci_hg38_chr.zip "https://ora.ox.ac.uk/objects/uuid:08e24957-7e76-438a-bd38-66c48008cf52/download_file?file_format=&safe_filename=1000G_loci_hg38_chr.zip&type_of_work=Dataset" && \
  	unzip 1000G_loci_hg38_chr.zip
  	rm 1000G_loci_hg38_chr.zip
  	sed -E -i 's|^X|chrX|' 1000G_loci_hg38/1kg.phase3.v5a_GRCh38nounref_loci_chrX.txt

  	# Probloci
  	wget -q -O probloci_chr.zip "https://ora.ox.ac.uk/objects/uuid:08e24957-7e76-438a-bd38-66c48008cf52/download_file?file_format=&safe_filename=probloci_chr.zip&type_of_work=Dataset" && \
  	unzip probloci_chr.zip
  	rm probloci_chr.zip

  	# Beagle5
  	wget -q -O beagle_chr.zip "https://ora.ox.ac.uk/objects/uuid:08e24957-7e76-438a-bd38-66c48008cf52/download_file?file_format=&safe_filename=beagle_chr.zip&type_of_work=Dataset" && \
  	unzip beagle_chr.zip
  	mv beagle/ beagle5/ && \
  	cd beagle5/ && \
  	mkdir -p tmp/ && \
  	mv chr*.1kg.phase3.v5a_GRCh38nounref.vcf.gz tmp/ && \
  	for i in {1..22} X; do zcat "tmp/chr${i}.1kg.phase3.v5a_GRCh38nounref.vcf.gz" | sed -E 's|^'${i}'|chr'${i}'|' | gzip > "chr${i}.1kg.phase3.v5a_GRCh38nounref.vcf.gz"; done && \
  	rm -rf tmp/ && \
  	cd ../
}

#Function call
battenbergRefDownload