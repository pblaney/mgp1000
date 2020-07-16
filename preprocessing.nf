// Myeloma Genome Project 1000
// Comprehensive pipeline for analysis of matched T/N Multiple Myeloma WGS data
// https://github.com/pblaney/mgp1000

// This portion of the pipeline is used for consistent preprocessing of all input WGS files.
// Both FASTQ and BAM files are supported formats for the input WGS files.
// The pipeline assumes that the input BAM reads were trimmed before alignment and that all FASTQs are in raw form.

log.info ''
log.info '##### Myeloma Genome Project 1000 Pipeline #####'
log.info '################################################'
log.info '~~~~~~~~~~~~~~~~~ PREPROCESSING ~~~~~~~~~~~~~~~~'
log.info '################################################'
log.info ''

def helpMessage() {
	log.info"""

	Usage:

		nextflow run preprocessing.nf -bg --input_format fastq --singularity_module "singularity/3.1" -profile preprocessing

	Mandatory Arguments:
		--input_format			       [str]  Format of input files, either: fastq or bam
		-profile				       [str]  Configuration profile to use, each profile described in nextflow.config file
									          Currently available: preprocessing

	Main Options:
		-bg 					      [flag]  Runs the pipeline processes in the background, this option should be included if deploying
		                                      pipeline with real data set so processes will not be cut if user disconnects from deployment
		                                      environment
		--singularity_module    [quoted str]  Indicates the name of the Singularity software module to be loaded for use in the pipeline,
		                                      this option is not needed if Singularity is natively installed on the deployment environment
		--skip_to_qc                   [str]  Skips directly to final Preprocessing QC step, either: yes or no 
		                                      can only be used in conjunction with bam as the input_format, should only be used for extreme
		                                      coverage BAMs that have been previously aligned with BWA MEM to the hg38 reference genome and
		                                      have adequate provenance to reflect this
		--help                        [flag]  Prints this message

	""".stripIndent()
}

// Print help message if requested
if( params.help ) exit 0, helpMessage()


// #################################################### \\
// ~~~~~~~~~~~~~ PARAMETER CONFIGURATION ~~~~~~~~~~~~~~ \\

// Declare the defaults for all pipeline parameters
params.input_format = "bam"
params.output_dir = "output"
params.skip_to_qc = "no"
params.help = null

// Set channels for reference files
Channel
	.fromPath( 'references/trimmomaticContaminants.fa' )
	.set{ trimmomatic_contaminants }

Channel
	.fromPath( 'references/hg38' )
	.set{ bwa_reference_dir }

Channel
	.fromPath( 'references/hg38/Homo_sapiens_assembly38.fasta' )
	.into{ reference_genome_fasta_forBaseRecalibrator;
	       reference_genome_fasta_forApplyBqsr;
	       reference_genome_fasta_forCollectWgsMetrics;
	       reference_genome_fasta_forCollectGcBiasMetrics }

Channel
	.fromPath( 'references/hg38/Homo_sapiens_assembly38.fasta.fai' )
	.into{ reference_genome_fasta_index_forBaseRecalibrator;
	       reference_genome_fasta_index_forApplyBqsr;
	       reference_genome_fasta_index_forCollectWgsMetrics;
	       reference_genome_fasta_index_forCollectGcBiasMetrics }

Channel
	.fromPath( 'references/hg38/Homo_sapiens_assembly38.dict' )
	.into{ reference_genome_fasta_dict_forBaseRecalibrator;
	       reference_genome_fasta_dict_forApplyBqsr;
	       reference_genome_fasta_dict_forCollectWgsMetrics;
	       reference_genome_fasta_dict_forCollectGcBiasMetrics }

Channel
	.fromPath( 'references/hg38/wgs_calling_regions.hg38.interval_list' )
	.set{ gatk_bundle_wgs_interval_list }

Channel
	.fromPath( 'references/hg38/wgs_calling_regions.hg38.bed' )
	.set{ gatk_bundle_wgs_bed}

Channel
	.fromPath( 'references/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz' )
	.set{ gatk_bundle_mills_1000G }

Channel
	.fromPath( 'references/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi' )
	.set{ gatk_bundle_mills_1000G_index }

Channel
	.fromPath( 'references/hg38/Homo_sapiens_assembly38.known_indels.vcf.gz' )
	.set{ gatk_bundle_known_indels }

Channel
	.fromPath( 'references/hg38/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi' )
	.set{ gatk_bundle_known_indels_index }

Channel
	.fromPath( 'references/hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz' )
	.set{ gatk_bundle_dbsnp138 }

Channel
	.fromPath( 'references/hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi' )
	.set{ gatk_bundle_dbsnp138_index }


// #################################################### \\
// ~~~~~~~~~~~~~~~~ PIPELINE PROCESSES ~~~~~~~~~~~~~~~~ \\

// Set the channel for all input BAM files
Channel
	.fromPath( 'input/*.bam' )
	.into{ input_mapped_bams;
	       input_mapped_bams_forQaulimap }

// GATK RevertSam ~ convert input mapped BAM files to unmapped BAM files
process revertMappedBam_gatk {
	publishDir "${params.output_dir}/preprocessing/revertedBams", mode: 'symlink'

	input:
	path bam_mapped from input_mapped_bams

	output:
	path bam_unmapped into unmapped_bams

	when:
	params.input_format == "bam"
	params.skip_to_qc == "no"

	script:
	bam_unmapped = "${bam_mapped}".replaceFirst(/\.bam/, ".unmapped.bam")
	"""
	gatk RevertSam \
	--java-options "-Xmx36G -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1 -XX:+AggressiveOpts" \
	--VERBOSITY ERROR \
	--MAX_RECORDS_IN_RAM 8000000 \
	-I "${bam_mapped}" \
	-O "${bam_unmapped}"
	"""
}

// biobambam bamtofastq ~ convert unmapped BAM files to paired FASTQ files
process bamToFastq_biobambam {
	publishDir "${params.output_dir}/preprocessing/biobambam", mode: 'symlink'

	input:
	path bam_unmapped from unmapped_bams

	output:
	tuple path(fastq_R1), path(fastq_R2) into converted_fastqs_forFastqc, converted_fastqs_forAlignment

	when:
	params.input_format == "bam"
	params.skip_to_qc == "no"

	script:
	fastq_R1 = "${bam_unmapped}".replaceFirst(/\.unmapped\.bam/, "_R1.fastq.gz")
	fastq_R2 = "${bam_unmapped}".replaceFirst(/\.unmapped\.bam/, "_R2.fastq.gz")
	"""
	bamtofastq \
	filename="${bam_unmapped}" \
	F="${fastq_R1}" F2="${fastq_R2}" \
	gz=1
	"""
}

// If input files are FASTQs, set channel up for both R1 and R2 reads then merge into single channel
input_R1_fastqs = Channel.fromPath( 'input/*R1*.f*q*' )
input_R2_fastqs = Channel.fromPath( 'input/*R2*.f*q*' )
input_fastqs = input_R1_fastqs.merge( input_R2_fastqs )

// Trimmomatic ~ trim low quality bases and clip adapters from reads
process fastqTrimming_trimmomatic {
	publishDir "${params.output_dir}/preprocessing/trimmomatic/trimmedFastqs", mode: 'symlink', pattern: "*${fastq_R1_trimmed}"	
	publishDir "${params.output_dir}/preprocessing/trimmomatic/trimmedFastqs", mode: 'symlink', pattern: "*${fastq_R2_trimmed}"
	publishDir "${params.output_dir}/preprocessing/trimmomatic/trimLogs", mode: 'move', pattern: "*${fastq_trim_log}"

	input:
	tuple path(input_R1_fastqs), path(input_R2_fastqs) from input_fastqs
	path trimmomatic_contaminants

	output:
	tuple path(fastq_R1_trimmed), path(fastq_R2_trimmed) into trimmed_fastqs_forFastqc, trimmed_fastqs_forAlignment
	path fastq_trim_log

	when:
	params.input_format == "fastq"

	script:
	fastq_R1_trimmed = "${input_R1_fastqs}".replaceFirst(/\.f.*q*/, ".trim.fastq.gz")
	fastq_R2_trimmed = "${input_R2_fastqs}".replaceFirst(/\.f.*q*/, ".trim.fastq.gz")
	fastq_R1_unpaired = "${input_R1_fastqs}".replaceFirst(/\.f.*q*/, ".unpaired.fastq.gz")
	fastq_R2_unpaired = "${input_R2_fastqs}".replaceFirst(/\.f.*q*/, ".unpaired.fastq.gz")
	fastq_trim_log = "${input_R1_fastqs}".replaceFirst(/[_\.]R1.*\.f.*q*/, ".trim.log")
	"""
	trimmomatic PE -threads 4 \
	"${input_R1_fastqs}" "${input_R2_fastqs}" \
	"${fastq_R1_trimmed}" "${fastq_R1_unpaired}" \
	"${fastq_R2_trimmed}" "${fastq_R2_unpaired}" \
	ILLUMINACLIP:${trimmomatic_contaminants}:2:30:10:1:true TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:35 \
	2> "${fastq_trim_log}"
	"""
}

// Depending on which input data type was used, set an input variable for the FastQC process
if( params.input_format == "bam" ) {
	fastqs_forFastqc = converted_fastqs_forFastqc
}
else {
	fastqs_forFastqc = trimmed_fastqs_forFastqc
}

// FastQC ~ generate sequence quality metrics for input FASTQ files
process fastqQualityControlMetrics_fastqc {
	publishDir "${params.output_dir}/preprocessing/fastqc", mode: 'move'

	input:
	tuple path(fastq_R1), path(fastq_R2) from fastqs_forFastqc

	output:
	tuple path(fastqc_R1_html), path(fastqc_R2_html)
	tuple path(fastqc_R1_zip), path(fastqc_R2_zip)

	when:
	params.skip_to_qc == "no"

	script:
	fastqc_R1_html = "${fastq_R1}".replaceFirst(/\.*fastq.gz/, "_fastqc.html")
	fastqc_R1_zip = "${fastq_R1}".replaceFirst(/\.*fastq.gz/, "_fastqc.zip")
	fastqc_R2_html = "${fastq_R2}".replaceFirst(/\.*fastq.gz/, "_fastqc.html")
	fastqc_R2_zip = "${fastq_R2}".replaceFirst(/\.*fastq.gz/, "_fastqc.zip")
	"""
	fastqc --outdir . "${fastq_R1}"
	fastqc --outdir . "${fastq_R2}"
	"""
}

// Depending on which input data type was used, set an input variable for the BWA process
if( params.input_format == "bam" ) {
	fastqs_forAlignment = converted_fastqs_forAlignment
}
else {
	fastqs_forAlignment = trimmed_fastqs_forAlignment
}

// BWA MEM + Sambamba ~ align trimmed FASTQ files to reference genome
process alignment_bwa {
	publishDir "${params.output_dir}/preprocessing/alignment/rawBams", mode: 'symlink', pattern: "*${bam_aligned}"
	publishDir "${params.output_dir}/preprocessing/alignment/flagstatLogs", mode: 'move', pattern: "*${bam_flagstat_log}"

	input:
	tuple path(fastq_R1), path(fastq_R2) from fastqs_forAlignment
	path bwa_reference_dir

	output:
	path bam_aligned into aligned_bams
	path bam_flagstat_log

	when:
	params.skip_to_qc == "no"

	script:
	sample_id = "${fastq_R1}".replaceFirst(/[_\.]R1.*\.f.*q*/, "")
	bam_aligned = "${sample_id}.bam"
	bam_flagstat_log = "${sample_id}.flagstat.log"
	"""
	bwa mem \
	-M -v 1 \
	-t 6 \
	-R '@RG\\tID:${sample_id}\\tSM:${sample_id}\\tLB:${sample_id}\\tPL:ILLUMINA' \
	"${bwa_reference_dir}/Homo_sapiens_assembly38.fasta" \
	"${fastq_R1}" "${fastq_R2}" \
	| \
	sambamba view \
	--sam-input \
	--nthreads=6 \
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

	when:
	params.skip_to_qc == "no"

	script:
	bam_fixed_mate = "${bam_aligned}".replaceFirst(/\.bam/, ".fixedmate.bam")
	"""
	gatk FixMateInformation \
	--java-options "-Xmx36G -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1 -XX:+AggressiveOpts" \
	--VERBOSITY ERROR \
	--VALIDATION_STRINGENCY SILENT \
	--ADD_MATE_CIGAR true \
	--MAX_RECORDS_IN_RAM 8000000 \
	--SORT_ORDER coordinate \
	-I "${bam_aligned}" \
	-O "${bam_fixed_mate}"
	"""
}

// Sambamba ~ mark duplicate alignments, remove them, and create BAM index
process markDuplicatesAndIndex_sambamba {
	publishDir "${params.output_dir}/preprocessing/markedDuplicates/bamsWithIndcies", mode: 'symlink'
	publishDir "${params.output_dir}/preprocessing/markedDuplicates/flagstatLogs", mode: 'move', pattern: "*${bam_markdup_flagstat_log}"

	input:
	path bam_fixed_mate from fixed_mate_bams

	output:
	path bam_marked_dup into marked_dup_bams_forDownsampleBam, marked_dup_bams_forApplyBqsr
	path bam_marked_dup_index
	path bam_markdup_flagstat_log

	when:
	params.skip_to_qc == "no"

	script:
	bam_marked_dup = "${bam_fixed_mate}".replaceFirst(/\.fixedmate\.bam/, ".markdup.bam")
	bam_marked_dup_index = "${bam_marked_dup}.bai"
	markdup_output_log = "${bam_fixed_mate}".replaceFirst(/\.fixedmate\.bam/, ".markdup.log")
	bam_markdup_flagstat_log = "${bam_fixed_mate}".replaceFirst(/\.fixedmate\.bam/, ".markdup.flagstat.log")
	"""
	sambamba markdup \
	--remove-duplicates \
	--nthreads 2 \
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

	when:
	params.skip_to_qc == "no"

	script:
	bam_marked_dup_downsampled = "${bam_marked_dup}".replaceFirst(/\.markdup\.bam/, ".markdup.downsampled.bam")
	"""
	gatk DownsampleSam \
	--java-options "-Xmx36G -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1 -XX:+AggressiveOpts" \
	--VERBOSITY ERROR \
	--MAX_RECORDS_IN_RAM 8000000 \
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
gatk_bundle_wgs_interval_list.combine( gatk_bundle_mills_1000G )
	.combine( gatk_bundle_mills_1000G_index )
	.combine( gatk_bundle_known_indels )
	.combine( gatk_bundle_known_indels_index )
	.combine( gatk_bundle_dbsnp138 )
	.combine( gatk_bundle_dbsnp138_index )
	.set{ gatk_reference_bundle }

reference_genome_fasta_forBaseRecalibrator.combine( reference_genome_fasta_index_forBaseRecalibrator )
	.combine( reference_genome_fasta_dict_forBaseRecalibrator )
	.set{ reference_genome_bundle_forBaseRecalibrator }

// GATK BaseRecalibrator ~ generate base quality score recalibration table based on covariates
process baseRecalibrator_gatk {
	publishDir "${params.output_dir}/preprocessing/baseRecalibration", mode: 'symlink'

	input:
	path bam_marked_dup_downsampled from downsampled_makred_dup_bams
	tuple path(reference_genome_fasta_forBaseRecalibrator), path(reference_genome_fasta_index_forBaseRecalibrator), path(reference_genome_fasta_dict_forBaseRecalibrator) from reference_genome_bundle_forBaseRecalibrator
	tuple path(gatk_bundle_wgs_interval_list), path(gatk_bundle_mills_1000G), path(gatk_bundle_mills_1000G_index), path(gatk_bundle_known_indels), path(gatk_bundle_known_indels_index), path(gatk_bundle_dbsnp138), path(gatk_bundle_dbsnp138_index) from gatk_reference_bundle

	output:
	path bqsr_table into base_quality_score_recalibration_data

	when:
	params.skip_to_qc == "no"

	script:
	bqsr_table = "${bam_marked_dup_downsampled}".replaceFirst(/\.markdup\.downsampled\.bam/, ".recaldata.table")
	"""
	gatk BaseRecalibrator \
	--java-options "-Xmx4G -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1 -XX:+AggressiveOpts" \
	--verbosity ERROR \
	--read-filter GoodCigarReadFilter \
	--reference "${reference_genome_fasta_forBaseRecalibrator}" \
	-L "${gatk_bundle_wgs_interval_list}" \
	-I "${bam_marked_dup_downsampled}" \
	-O "${bqsr_table}" \
	--known-sites "${gatk_bundle_mills_1000G}" \
	--known-sites "${gatk_bundle_known_indels}" \
	--known-sites "${gatk_bundle_dbsnp138}"
	"""
}

// Create additional channel for the reference FASTA to be used in GATK ApplyBQSR process
reference_genome_fasta_forApplyBqsr.combine( reference_genome_fasta_index_forApplyBqsr )
	.combine( reference_genome_fasta_dict_forApplyBqsr )
	.set{ reference_genome_bundle_forApplyBqsr }

// GATK ApplyBQSR ~ apply base quality score recalibration using generated table
process applyBqsr_gatk {
	publishDir "${params.output_dir}/preprocessing/finalPreprocessedBams", mode: 'copy', pattern: '*.{bam,bai}'

	input:
	path bam_marked_dup from marked_dup_bams_forApplyBqsr
	tuple path(reference_genome_fasta_forApplyBqsr), path(reference_genome_fasta_index_forApplyBqsr), path(reference_genome_fasta_dict_forApplyBqsr) from reference_genome_bundle_forApplyBqsr
	path bqsr_table from base_quality_score_recalibration_data

	output:
	path bam_preprocessed_final into final_preprocessed_bams_forCollectWgsMetrics, final_preprocessed_bams_forCollectGcBiasMetrics

	when:
	params.skip_to_qc == "no"

	script:
	bam_preprocessed_final = "${bam_marked_dup}".replaceFirst(/\.markdup\.bam/, ".final.bam")
	"""
	gatk ApplyBQSR \
	--java-options "-Xmx4G -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1 -XX:+AggressiveOpts" \
	--verbosity ERROR \
	--read-filter GoodCigarReadFilter \
	--reference "${reference_genome_fasta_forApplyBqsr}" \
	-I "${bam_marked_dup}" \
	-O "${bam_preprocessed_final}" \
	--bqsr-recal-file "${bqsr_table}"
	"""
}

// Create additional channel for the reference FASTA and interfal list to be used in GATK CollectWgsMetrics process
reference_genome_fasta_forCollectWgsMetrics.combine( reference_genome_fasta_index_forCollectWgsMetrics )
	.combine( reference_genome_fasta_dict_forCollectWgsMetrics )
	.set{ reference_genome_bundle_forCollectWgsMetrics }

// GATK CollectWgsMetrics ~ generate covearge and performance metrics from final BAM
process collectWgsMetrics_gatk {
	publishDir "${params.output_dir}/preprocessing/coverageMetrics", mode: 'move'

	input:
	path bam_preprocessed_final from final_preprocessed_bams_forCollectWgsMetrics
	tuple path(reference_genome_fasta_forCollectWgsMetrics), path(reference_genome_fasta_index_forCollectWgsMetrics), path(reference_genome_fasta_dict_forCollectWgsMetrics) from reference_genome_bundle_forCollectWgsMetrics

	output:
	path coverage_metrics

	when:
	params.skip_to_qc == "no"

	script:
	sample_id = "${bam_preprocessed_final}".replaceFirst(/\.final\.bam/, "")
	coverage_metrics = "${sample_id}.coverage.metrics.txt"
	"""
	gatk CollectWgsMetrics \
	--java-options "-Xmx6G -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1 -XX:+AggressiveOpts" \
	--VERBOSITY ERROR \
	--INCLUDE_BQ_HISTOGRAM \
	--MINIMUM_BASE_QUALITY 20 \
	--MINIMUM_MAPPING_QUALITY 20 \
	--REFERENCE_SEQUENCE "${reference_genome_fasta_forCollectWgsMetrics}" \
	-I "${bam_preprocessed_final}" \
	-O "${coverage_metrics}"
	"""
}

// Create additional channel for the reference FASTA and interfal list to be used in GATK CollectWgsMetrics process
reference_genome_fasta_forCollectGcBiasMetrics.combine( reference_genome_fasta_index_forCollectGcBiasMetrics )
	.combine( reference_genome_fasta_dict_forCollectGcBiasMetrics )
	.set{ reference_genome_bundle_forCollectGcBiasMetrics }

// GATK CollectGcBiasMetrics ~ generate GC content bias in reads in final BAM
process collectGcBiasMetrics_gatk {
	publishDir "${params.output_dir}/preprocessing/gcBiasMetrics", mode: 'move'

	input:
	path bam_preprocessed_final from final_preprocessed_bams_forCollectGcBiasMetrics
	tuple path(reference_genome_fasta_forCollectGcBiasMetrics), path(reference_genome_fasta_index_forCollectGcBiasMetrics), path(reference_genome_fasta_dict_forCollectGcBiasMetrics) from reference_genome_bundle_forCollectGcBiasMetrics

	output:
	path gc_bias_metrics
	path gc_bias_chart
	path gc_bias_summary

	when:
	params.skip_to_qc == "no"

	script:
	sample_id = "${bam_preprocessed_final}".replaceFirst(/\.final\.bam/, "")
	gc_bias_metrics = "${sample_id}.gcbias.metrics.txt"
	gc_bias_chart = "${sample_id}.gcbias.metrics.pdf"
	gc_bias_summary = "${sample_id}.gcbias.summary.txt"
	"""
	gatk CollectGcBiasMetrics \
	--java-options "-Xmx6G -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1 -XX:+AggressiveOpts" \
	--VERBOSITY ERROR \
	--REFERENCE_SEQUENCE "${reference_genome_fasta_forCollectGcBiasMetrics}" \
	-I "${bam_preprocessed_final}" \
	-O "${gc_bias_metrics}" \
	--CHART_OUTPUT "${gc_bias_chart}" \
	--SUMMARY_OUTPUT "${gc_bias_summary}"
	"""
}

// If input BAM files have excessively high coverage (>100x), were previously aligned with same reference genome,
// and aligner then there is the option to skip directly to QC processes

// Qualimap BAM QC ~ perform quality control analysis on extreme coverage (>100x) BAMs
process extremeBamQualityControl_qualimap {
	publishDir "${params.output_dir}/preprocessing/qualimapBamQc", mode: 'move'

	input:
	path bam_mapped_high_coverage from input_mapped_bams_forQaulimap

	output:
	path "${bam_mapped_high_coverage.baseName}"

	when:
	params.input_format == "bam"
	params.skip_to_qc == "yes"

	script:
	"""
	qualimap --java-mem-size=12G \
	bamqc \
	-bam "${bam_mapped_high_coverage}" \
	--paint-chromosome-limits \
	--genome-gc-distr HUMAN \
	-nt 8 \
	-outformat HTML \
	-outdir "${bam_mapped_high_coverage.baseName}"
	"""
}
