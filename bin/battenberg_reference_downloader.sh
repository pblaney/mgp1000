#!/usr/bin/env bash
# Download all the necessary files for use in Battenberg

# Function
battenbergRefDownload() {

	# GC Correction
	echo "Beginning GC Correction Download....."
	wget -q -O GC_correction_hg38_chr.zip "https://www.dropbox.com/sh/bize1n830t0mgzb/AADQD4DTJOF75YmhBDDoQ9nla/GC_correction_hg38?dl=0&lst=" && \
  	unzip GC_correction_hg38_chr.zip -x /
  	rm GC_correction_hg38_chr.zip && \
	mv 1000G_GC_chr*.txt.gz GC_correction_hg38/

	echo 
	echo "GC Correction  }}}---->>>  D O N E"
	echo 
	sleep 5

	# RT Correction
	echo "Beginning RT Correction Download....."
  	wget -q -O RT_correction_hg38.zip "https://www.dropbox.com/sh/bize1n830t0mgzb/AABZ2uM13YMYB_q6X1pP1McJa/RT_correction_hg38?dl=0&lst=" && \
  	unzip RT_correction_hg38.zip -x /
  	rm RT_correction_hg38.zip && \
  	mv 1000G_RT_chr*.txt.gz RT_correction_hg38/

  	echo 
  	echo "RT Correction  }}}---->>>  D O N E"
  	echo 
  	sleep 5
  	
  	# Shapeit2
  	echo "Beginning Shapeit2 Download....."
  	wget -q -O shapeit2_chr.zip "https://www.dropbox.com/sh/bize1n830t0mgzb/AAAbmJCGL04zOqNm7v7HRfjza/shapeit2?dl=0&lst=" && \
  	unzip shapeit2_chr.zip -x /
  	rm shapeit2_chr.zip && \
  	mv ALL.v1a.shapeit2_integrated_chr*.GRCh38.20181129.phased.hap.gz shapeit2/ && \
  	mv ALL.v1a.shapeit2_integrated_chr*.GRCh38.20181129.phased.legend temp/
  	for i in {1..22}; do cat "temp/ALL.v1a.shapeit2_integrated_chr${i}.GRCh38.20181129.phased.legend" | sed -E 's|^'${i}'|chr'${i}'|' > "ALL.v1a.shapeit2_integrated_chr${i}.GRCh38.20181129.phased.legend"; done
  	for i in X_PAR1 X_PAR2 X_nonPAR; do cat "temp/ALL.v1a.shapeit2_integrated_chr${i}.GRCh38.20181129.phased.legend" | sed -E 's|^X|chrX|' > "ALL.v1a.shapeit2_integrated_chr${i}.GRCh38.20181129.phased.legend"; done
  	mv ALL.v1a.shapeit2_integrated_chr*.GRCh38.20181129.phased.legend shapeit2/
  	rm temp/*

  	echo 
  	echo "Shapeit2  }}}---->>>  D O N E"
  	echo 
  	sleep 5

  	# Impute
  	echo "Beginning Imputation Download....."
  	wget -q -O imputation_chr.zip "https://www.dropbox.com/sh/bize1n830t0mgzb/AABkSpOJdQthXn5ELVZatnrxa/imputation?dl=0&lst=" && \
	unzip imputation_chr.zip -x /
	rm imputation_chr.zip && \
	mv genetic_map_chr*_combined_b38.txt imputation/
	mv impute_info.txt temp/ && \
	awk 'BEGIN {OFS="\t"} {print "chr"$1,$2,$3,$4,$5,$6}' temp/impute_info.txt > impute_info.txt
	rm temp/*

	echo 
	echo "Imputation  }}}---->>>  D O N E"
	echo 
	sleep 5

  	# 1000G
  	echo "Beginning 1000 Genome Loci Download....."
  	wget -q -O 1000G_loci_hg38_chr.zip "https://www.dropbox.com/sh/bize1n830t0mgzb/AAD8szWaYFjkeFElpIn9Kxcra/1000G_loci_hg38?dl=0&lst=" && \
  	unzip 1000G_loci_hg38_chr.zip -x /
  	rm 1000G_loci_hg38_chr.zip && \
  	rm 1kg.phase3.v5a_GRCh38nounref_loci_chrstring_chr*.txt && \
  	mv 1kg.phase3.v5a_GRCh38nounref_allele_index_chr*.txt 1000G_loci_hg38/ && \
  	mv 1kg.phase3.v5a_GRCh38nounref_loci_chr*.txt temp/
  	for i in {1..22} X; do cat "temp/1kg.phase3.v5a_GRCh38nounref_loci_chr${i}.txt" | sed -E 's|^'${i}'|chr'${i}'|' > "1kg.phase3.v5a_GRCh38nounref_loci_chr${i}.txt"; done
	mv 1kg.phase3.v5a_GRCh38nounref_loci_chr*.txt 1000G_loci_hg38/ && \
	rm temp/*

	echo 
	echo "1000 Genome Loci  }}}---->>>  D O N E"
	echo 
	sleep 5

  	# Probloci
  	echo "Beginning Problem Loci Download....."
  	wget -q -O probloci_chr.zip "https://www.dropbox.com/sh/bize1n830t0mgzb/AAByAZhtoyceFWiLHQ4VVVHLa/probloci?dl=0&lst=" && \
  	unzip probloci_chr.zip -x /
  	rm probloci_chr.zip && \
  	mv probloci.txt.gz probloci/

  	echo 
  	echo "Problem Loci  }}}---->>>  D O N E"
  	echo 
  	sleep 5

  	# Beagle5
  	echo "Beginning Beagle5 Download....."
  	wget -q -O beagle_chr.zip "https://www.dropbox.com/sh/bize1n830t0mgzb/AADTvENGArNqA72vUf_OHRhza/beagle5?dl=0&lst=" && \
  	unzip beagle_chr.zip -x /
  	rm beagle_chr.zip && \
  	rm combine_chrX_plink.R && \
  	mv chr*.1kg.phase3.v5a_GRCh38nounref.vcf.gz temp/
  	for i in {1..22} X; do zcat "temp/chr${i}.1kg.phase3.v5a_GRCh38nounref.vcf.gz" | sed -E 's|^'${i}'|chr'${i}'|' | gzip > "chr${i}.1kg.phase3.v5a_GRCh38nounref.vcf.gz"; done
  	mv chr*.1kg.phase3.v5a_GRCh38nounref.vcf.gz beagle5/ && \
  	rm -rf temp/
  	mv plink.chr*.GRCh38.map temp/ && \
  	for i in {1..22} X; do zcat "temp/plink.chr${i}.GRCh38.map" | sed -E 's|^'${i}'|chr'${i}'|' > "plink.chr${i}.GRCh38.map"; done
  	mv plink.chr*.GRCh38.map beagle5/ && \
  	rm -rf temp/

  	echo 
  	echo "Beagle5  }}}---->>>  D O N E"
  	echo 
  	sleep 5

}

#Function call
battenbergRefDownload

echo "D O N E"
echo 