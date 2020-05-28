// ################################################################################
// |\       /|   /---------\  /---------\    /|   /------\    /------\    /------\
// | \     / |  |             |         |   / |  |        |  |        |  |        |
// |  \   /  |  |     |----\  |---------/     |  |        |  |        |  |        |
// |   \ /   |  |          |  |               |  |        |  |        |  |        |
// |         |   \---------/  |               |   \------/    \------/    \------/
// ################################################################################

// Myeloma Genome Project 1000 
// Comprehensive pipeline for analysis of matched T/N Multiple Myeloma WGS data
// https://github.com/pblaney/mgp1000


// ################################################ \\
//                                                  \\
// ~~~~~~~~~~~~~~~~ CONFIGURATION ~~~~~~~~~~~~~~~~~ \\
//                                                  \\
// ################################################ \\


// NEED TO CHANGE MODE TO COPY WHEN FINALIZED

// Declare the defaults for all pipeline parameters
params.input_format = "bam"
params.output_dir = "output"

// Set path to reference files
trimmomatic_contaminants = Channel.fromPath( 'references/trimmomatic.fa' )
bwa_reference_dir = Channel.fromPath( 'references/hg38/bwa' )
reference_genome_fasta = Channel.fromPath( 'references/hg38/bwa/genome.fa' )
reference_genome_fasta_index = Channel.fromPath( 'references/hg38/bwa/genome.fa.fai' )
reference_genome_fasta_dict = Channel.fromPath( 'references/hg38/bwa/genome.dict' )
gatk_bundle_wgs_interval_list = Channel.fromPath( 'references/hg38/gatkBundle/wgs_calling_regions.hg38.interval_list' )
gatk_bundle_mills_1000G = Channel.fromPath( 'references/hg38/gatkBundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz' )
gatk_bundle_mills_1000G_index = Channel.fromPath( 'references/hg38/gatkBundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi' )
gatk_bundle_known_indels = Channel.fromPath( 'references/hg38/gatkBundle/Homo_sapiens_assembly38.known_indels.vcf.gz' )
gatk_bundle_known_indels_index = Channel.fromPath( 'references/hg38/gatkBundle/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi' )
gatk_bundle_dbsnp138 = Channel.fromPath( 'references/hg38/gatkBundle/Homo_sapiens_assembly38.dbsnp138.vcf' )
gatk_bundle_dbsnp138_index = Channel.fromPath( 'references/hg38/gatkBundle/Homo_sapiens_assembly38.dbsnp138.vcf.idx' )


// ################################################ \\
//                                                  \\
// ~~~~~~~~~~~~~~~~ PRE-PROCESSING ~~~~~~~~~~~~~~~~ \\
//                                                  \\
// ################################################ \\

// Set the path to all input BAM files
input_bams = Channel.fromPath( 'input/*.bam' )

// biobambam ~ convert any input BAM files to FASTQ files
process bamToFastq_biobambam {
	publishDir "${params.output_dir}/preprocessing/biobambam", mode: 'symlink'

	input:
	path bam from input_bams

	output:
	tuple path(fastq_R1), path(fastq_R2) into converted_fastqs

	when:
	params.input_format == "bam"

	script:
	fastq_R1 = "${bam}".replaceFirst(/.bam$/, "_R1.fastq.gz")
	fastq_R2 = "${bam}".replaceFirst(/.bam$/, "_R2.fastq.gz")
	"""
	bamtofastq filename="${bam}" F="${fastq_R1}" F2="${fastq_R2}" gz=1
	"""
}

// Set the path to all input FASTQ files. If started with BAM files, use the converted FASTQs
if( params.input_format == "bam" ) {
	input_fastqs = converted_fastqs
}
else {
	input_R1_fastqs = Channel.fromPath( 'input/*_R1.f*q*' )
	input_R2_fastqs = Channel.fromPath( 'input/*_R2.f*q*' )
	input_fastqs = input_R1_fastqs.merge( input_R2_fastqs )
}

// Trimmomatic ~ trim low quality bases and clip adapters from reads
process fastqTrimming_trimmomatic {
	publishDir "${params.output_dir}/preprocessing/trimmomatic/trimmedFastqs", mode: 'symlink', pattern: "*${fastq_R1_trimmed}"	
	publishDir "${params.output_dir}/preprocessing/trimmomatic/trimmedFastqs", mode: 'symlink', pattern: "*${fastq_R2_trimmed}"
	publishDir "${params.output_dir}/preprocessing/trimmomatic/trimLogs", mode: 'symlink', pattern: "*${fastq_trim_log}"

	input:
	tuple path(input_R1_fastqs), path(input_R2_fastqs) from input_fastqs
	path trimmomatic_contaminants

	output:
	tuple path(fastq_R1_trimmed), path(fastq_R2_trimmed) into trimmed_fastqs_forFastqc, trimmed_fastqs_forAlignment
	path fastq_trim_log

	script:
	fastq_R1_trimmed = "${input_R1_fastqs}".replaceFirst(/.fastq.gz$/, ".trim.fastq.gz")
	fastq_R2_trimmed = "${input_R2_fastqs}".replaceFirst(/.fastq.gz$/, ".trim.fastq.gz")
	fastq_R1_unpaired = "${input_R1_fastqs}".replaceFirst(/.fastq.gz$/, ".unpaired.fastq.gz")
	fastq_R2_unpaired = "${input_R2_fastqs}".replaceFirst(/.fastq.gz$/, ".unpaired.fastq.gz")
	fastq_trim_log = "${input_R1_fastqs}".replaceFirst(/_R1.fastq.gz$/, ".trim.log")
	"""
	trimmomatic PE -threads 8 \
	"${input_R1_fastqs}" "${input_R2_fastqs}" \
	"${fastq_R1_trimmed}" "${fastq_R1_unpaired}" \
	"${fastq_R2_trimmed}" "${fastq_R2_unpaired}" \
	ILLUMINACLIP:${trimmomatic_contaminants}:2:30:10:1:true TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:35 \
	2> "${fastq_trim_log}"
	"""
}

// FastQC ~ generate sequence quality metrics for input FASTQ files
process fastqQualityControlMetrics_fastqc {
	publishDir "${params.output_dir}/preprocessing/fastqc", mode: 'symlink'

	input:
	tuple path(fastq_R1), path(fastq_R2) from trimmed_fastqs_forFastqc

	output:
	tuple path(fastqc_R1_html), path(fastqc_R2_html) into fastqc_reports
	tuple path(fastqc_R1_zip), path(fastqc_R2_zip)

	script:
	fastqc_R1_html = "${fastq_R1}".replaceFirst(/.fastq.gz$/, "_fastqc.html")
	fastqc_R1_zip = "${fastq_R1}".replaceFirst(/.fastq.gz$/, "_fastqc.zip")
	fastqc_R2_html = "${fastq_R2}".replaceFirst(/.fastq.gz$/, "_fastqc.html")
	fastqc_R2_zip = "${fastq_R2}".replaceFirst(/.fastq.gz$/, "_fastqc.zip")
	"""
	fastqc --outdir . "${fastq_R1}"
	fastqc --outdir . "${fastq_R2}"
	"""
}

// BWA MEM + Sambamba ~ align trimmed FASTQ files to reference genome
process alignment_bwa {
	publishDir "${params.output_dir}/preprocessing/alignment/rawBams", mode: 'symlink', pattern: "*${bam_aligned}"
	publishDir "${params.output_dir}/preprocessing/alignment/flagstatLogs", mode: 'symlink', pattern: "*${bam_flagstat_log}"

	input:
	tuple path(fastq_R1_trimmed), path(fastq_R2_trimmed) from trimmed_fastqs_forAlignment
	path bwa_reference_dir

	output:
	path bam_aligned into aligned_bams
	path bam_flagstat_log

	script:
	sample_id = "${fastq_R1_trimmed}".replaceFirst(/_R1.trim.fastq.gz$/, "")
	bam_aligned = "${sample_id}.bam"
	bam_flagstat_log = "${sample_id}.flagstat.log"
	"""
	bwa mem \
	-M -v 1 \
	-t 8 \
	-R '@RG\\tID:${sample_id}\\tSM:${sample_id}\\tLB:${sample_id}\\tPL:ILLUMINA' \
	"${bwa_reference_dir}/genome.fa" \
	"${fastq_R1_trimmed}" "${fastq_R2_trimmed}" \
	| \
	sambamba view \
	--sam-input \
	--nthreads=8 \
	--filter='mapping_quality>=10' \
	--format=bam \
	--compression-level=0 \
	--output-filename "${bam_aligned}" \
	/dev/stdin

	sambamba flagstat "${bam_aligned}" > "${bam_flagstat_log}"
	"""
}

// GATK FixMateInformation ~ veryify/fix mate-pair information and sort output BAM by coordinate
process fixMateInformationAndSort_gatk {
	publishDir "${params.output_dir}/preprocessing/fixMateBams", mode: 'symlink'

	input:
	path bam_aligned from aligned_bams

	output:
	path bam_fixed_mate into fixed_mate_bams

	script:
	bam_fixed_mate = "${bam_aligned}".replaceFirst(/.bam$/, ".fixedmate.bam")
	"""
	gatk FixMateInformation \
	--java-options "-Xmx24576m -XX:ParallelGCThreads=1" \
	--VALIDATION_STRINGENCY SILENT \
	--ADD_MATE_CIGAR true \
	--SORT_ORDER coordinate \
	-I "${bam_aligned}" \
	-O "${bam_fixed_mate}"
	"""
}

// Sambamba ~ mark duplicate alignments and create BAM index
process markDuplicatesAndIndex_sambamba {
	publishDir "${params.output_dir}/preprocessing/markedDuplicates/bamsWithIndcies", mode: 'symlink'
	publishDir "${params.output_dir}/preprocessing/markedDuplicates/flagstatLogs", mode: 'symlink', pattern: "*${bam_markdup_flagstat_log}"

	input:
	path bam_fixed_mate from fixed_mate_bams

	output:
	path bam_marked_dup into marked_dup_bams_forDownsampleBam, marked_dup_bams_forApplyBqsr
	path bam_marked_dup_index
	path bam_markdup_flagstat_log
	//into marked_dup_bam_indices

	script:
	bam_marked_dup = "${bam_fixed_mate}".replaceFirst(/.fixedmate.bam$/, ".markdup.bam")
	bam_marked_dup_index = "${bam_marked_dup}.bai"
	markdup_output_log = "${bam_fixed_mate}".replaceFirst(/.fixedmate.bam$/, ".markdup.log")
	bam_markdup_flagstat_log = "${bam_fixed_mate}".replaceFirst(/.fixedmate.bam$/, ".markdup.flagstat.log")
	"""
	sambamba markdup \
	--remove-duplicates \
	--nthreads 4 \
	--hash-table-size 525000 \
	--overflow-list-size 525000 \
	"${bam_fixed_mate}" \
	"${bam_marked_dup}" \
	2> "${markdup_output_log}"

	sambamba flagstat "${bam_marked_dup}" > "${bam_markdup_flagstat_log}"

	sambamba index "${bam_marked_dup}" "${bam_marked_dup_index}"
	"""	
}

// GATK DownsampleSam ~ downsample BAM file to use random subset for generating BSQR table
process downsampleBam_gatk {
	publishDir  "${params.output_dir}/preprocessing/downsampleBams", mode: 'symlink'

	input:
	path bam_marked_dup from marked_dup_bams_forDownsampleBam

	output:
	path bam_marked_dup_downsampled into downsampled_makred_dup_bams

	script:
	bam_marked_dup_downsampled = "${bam_marked_dup}".replaceFirst(/.markdup.bam$/, ".markdup.downsampled.bam")
	"""
	gatk DownsampleSam \
	--java-options "-Xmx80g -XX:ParallelGCThreads=1" \
	--STRATEGY Chained \
	--RANDOM_SEED 1000 \
	--CREATE_INDEX \
	--VALIDATION_STRINGENCY SILENT \
	-P 0.1 \
	-I "${bam_marked_dup}" \
	-O "${bam_marked_dup_downsampled}"
	"""
}

// Combine all needed GATK bundle files and reference FASTA into one channel for use in GATK BaseRecalibrator process
gatk_reference_bundle = gatk_bundle_wgs_interval_list.combine(gatk_bundle_mills_1000G)
							.combine(gatk_bundle_mills_1000G_index)
							.combine(gatk_bundle_known_indels)
							.combine(gatk_bundle_known_indels_index)
							.combine(gatk_bundle_dbsnp138)
							.combine(gatk_bundle_dbsnp138_index)
reference_genome_bundle = reference_genome_fasta.combine(reference_genome_fasta_index)
							.combine(reference_genome_fasta_dict)

// GATK BaseRecalibrator ~ generate base quality score recalibration table based on covariates
process baseRecalibrator_gatk {
	publishDir "${params.output_dir}/preprocessing/baseRecalibration", mode: 'symlink'

	input:
	path bam_marked_dup_downsampled from downsampled_makred_dup_bams
	tuple path(reference_genome_fasta), path(reference_genome_fasta_index), path(reference_genome_fasta_dict) from reference_genome_bundle
	tuple path(gatk_bundle_wgs_interval_list), path(gatk_bundle_mills_1000G), path(gatk_bundle_mills_1000G_index), path(gatk_bundle_known_indels), path(gatk_bundle_known_indels_index), path(gatk_bundle_dbsnp138), path(gatk_bundle_dbsnp138_index) from gatk_reference_bundle

	output:
	path bqsr_table into base_quality_score_recalibration_data

	script:
	bqsr_table = "${bam_marked_dup_downsampled}".replaceFirst(/.markdup.downsampled.bam$/, ".recaldata.table")
	"""
	gatk BaseRecalibrator \
	--java-options "-Xmx24576m -XX:ParallelGCThreads=1" \
	--verbosity ERROR \
	--reference "${reference_genome_fasta}" \
	-L "${gatk_bundle_wgs_interval_list}" \
	-I "${bam_marked_dup_downsampled}" \
	-O "${bqsr_table}" \
	--known-sites "${gatk_bundle_mills_1000G}" \
	--known-sites "${gatk_bundle_known_indels}" \
	--known-sites "${gatk_bundle_dbsnp138}"
	"""
}

// GATK ApplyBQSR ~ apply base quality score recalibration using generated table
process applyBqsr_gatk {
	publishDir "${params.output_dir}/preprocessing/finalPreprocessedBam", mode: 'symlink'

	input:
	path bam_marked_dup from marked_dup_bams_forApplyBqsr
	tuple path(reference_genome_fasta), path(reference_genome_fasta_index), path(reference_genome_fasta_dict) from reference_genome_bundle
	path bqsr_table from base_quality_score_recalibration_data

	output:
	path bam_preprocessed_final into final_preprocessed_bams

	script:
	bam_preprocessed_final = "${bam_marked_dup}".replaceFirst(/.markdup.bam$/, ".final.bam")
	"""
	gatk ApplyBQSR \
	--java-options "-Xmx24576m -XX:ParallelGCThreads=1" \
	--reference "${reference_genome_fasta}" \
	-I "${bam_marked_dup}" \
	-O "${bam_preprocessed_final}" \
	--bqsr-recal-file "${bqsr_table}"
	"""
}