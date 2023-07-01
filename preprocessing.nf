// Myeloma Genome Pipeline 1000
// Comprehensive pipeline for analysis of matched T/N Multiple Myeloma WGS data
// https://github.com/pblaney/mgp1000

// This module of the pipeline is used for consistent preprocessing of all input WGS files.
// Both FASTQ and BAM files are supported formats for the input WGS files.

import java.text.SimpleDateFormat;
def workflowTimestamp = "${workflow.start.format('MM-dd-yyyy HH:mm')}"

def helpMessage() {
	log.info"""
	                             .------------------------.
	                            |    .-..-. .--. .---.     |
	                            |    : `' :: .--': .; :    |
	                            |    : .. :: : _ :  _.'    |
	                            |    : :; :: :; :: :       |
	                            |    :_;:_;`.__.':_;       |
	                            |   ,-. .--.  .--.  .--.   |
	                            | .'  :: ,. :: ,. :: ,. :  |
	                            |   : :: :: :: :: :: :: :  |
	                            |   : :: :; :: :; :: :; :  |
	                            |   :_;`.__.'`.__.'`.__.'  |
	                             .________________________.

                  ░█▀█░█▀▄░█▀▀░█▀█░█▀▄░█▀█░█▀▀░█▀▀░█▀▀░█▀▀░▀█▀░█▀█░█▀▀
                  ░█▀▀░█▀▄░█▀▀░█▀▀░█▀▄░█░█░█░░░█▀▀░▀▀█░▀▀█░░█░░█░█░█░█
                  ░▀░░░▀░▀░▀▀▀░▀░░░▀░▀░▀▀▀░▀▀▀░▀▀▀░▀▀▀░▀▀▀░▀▀▀░▀░▀░▀▀▀

	Usage:
	  nextflow run preprocessing.nf --run_id STR --input_format STR -profile preprocessing
	  [-bg] [-resume] [--lane_split STR] [--input_dir PATH] [--output_dir PATH] [--email STR]
	  [--seq_protocol STR][--cpus INT] [--memory STR] [--queue_size INT] [--executor STR] [--help]

	Mandatory Arguments:
	  --run_id                       STR  Unique identifier for pipeline run
	  --input_format                 STR  Format of input files
	                                      [Default: fastq | Available: fastq, bam]
	  -profile                       STR  Configuration profile to use, must use preprocessing                               

	Main Options:
	  -bg                           FLAG  Runs the pipeline processes in the background, this
	                                      option should be included if deploying pipeline with
	                                      real data set so processes will not be cut if user
	                                      disconnects from deployment environment
	  -resume                       FLAG  Successfully completed tasks are cached so that if
	                                      the pipeline stops prematurely the previously
	                                      completed tasks are skipped while maintaining their
	                                      output
	  --lane_split                   STR  Determines if input FASTQs are lane split per R1/R2
	                                      [Default: no | Available: yes, no]
	  --input_dir                   PATH  Directory that holds BAMs and associated index files,
	                                      this should be given as an absolute path
	                                      [Default: input/]
	  --output_dir                  PATH  Directory that will hold all output files this should
	                                      be given as an absolute path
	                                      [Default: output/]
	  --email                        STR  Email address to send workflow completion/stoppage
	                                      notification
	  --seq_protocol                 STR  Sequencing protocol of the input, WGS for whole-genome
	                                      and WES for whole-exome
	                                      [Default: WGS | Available: WGS, WES]
	  --cpus                         INT  Globally set the number of cpus to be allocated
	  --memory                       STR  Globally set the amount of memory to be allocated,
	                                      written as '##.GB' or '##.MB'
	  --queue_size                   INT  Set max number of tasks the pipeline will launch
	                                      [Default: 100]
	  --executor                     STR  Set the job executor for the run
	                                      [Default: slurm | Available: local, slurm, lsf]
	  --help                        FLAG  Prints this message
	""".stripIndent()
}

// #################################################### \\
// ~~~~~~~~~~~~~ PARAMETER CONFIGURATION ~~~~~~~~~~~~~~ \\

// Declare the defaults for all pipeline parameters
params.input_dir = "${workflow.projectDir}/input"
params.output_dir = "${workflow.projectDir}/output"
params.run_id = null
params.seq_protocol = "WGS"
params.input_format = "fastq"
params.lane_split = "no"
params.email = null
params.trimmomatic_min_len = 35
params.skip_trimming = "no"
params.cpus = null
params.memory = null
params.queue_size = 100
params.executor = 'slurm'
params.help = null

// Print help message if requested
if( params.help ) exit 0, helpMessage()

// Print erro message if user-defined input/output directories does not exist
if( !file(params.input_dir).exists() ) exit 1, "The user-specified input directory does not exist in filesystem."

// Print error messages if required parameters are not set
if( params.run_id == null ) exit 1, "The run command issued does not have the '--run_id' parameter set. Please set the '--run_id' parameter to a unique identifier for the run."

if( params.input_format == null ) exit 1, "The run command issued does not have the '--input_format' parameter set. Please set the '--input_format' parameter to either bam or fastq depending on input data."

if( params.input_format == "bam" & params.skip_trimming == "yes" ) exit 1, "This run command cannot be executed. If '--input_format' parameter is 'bam', then trimming must be performed. Please set '--skip_trimming to 'no'."

// Set channels for reference files
Channel
	.fromPath( 'references/trimmomaticContaminants.fa' )
	.set{ trimmomatic_contaminants }

Channel
	.fromPath( 'references/hg38' )
	.set{ bwa_reference_dir }

Channel
	.fromPath( 'references/hg38/Homo_sapiens_assembly38.fasta' )
	.into{ reference_genome_fasta_forRealignment;
		   reference_genome_fasta_forBaseRecalibrator;
	       reference_genome_fasta_forApplyBqsr;
	       reference_genome_fasta_forCollectWgsMetrics;
	       reference_genome_fasta_forCollectGcBiasMetrics;
	       reference_genome_fasta_forCollectHsMetrics }

Channel
	.fromPath( 'references/hg38/Homo_sapiens_assembly38.fasta.fai' )
	.into{ reference_genome_fasta_index_forRealignment;
	       reference_genome_fasta_index_forBaseRecalibrator;
	       reference_genome_fasta_index_forApplyBqsr;
	       reference_genome_fasta_index_forCollectWgsMetrics;
	       reference_genome_fasta_index_forCollectGcBiasMetrics;
	       reference_genome_fasta_index_forCollectHsMetrics }

Channel
	.fromPath( 'references/hg38/Homo_sapiens_assembly38.dict' )
	.into{ reference_genome_fasta_dict_forRealignment;
	       reference_genome_fasta_dict_forBaseRecalibrator;
	       reference_genome_fasta_dict_forApplyBqsr;
	       reference_genome_fasta_dict_forCollectWgsMetrics;
	       reference_genome_fasta_dict_forCollectGcBiasMetrics;
	       reference_genome_fasta_dict_forCollectHsMetrics }

if( params.seq_protocol == "WGS" ) {
    Channel
        .value( file('references/hg38/wgs_calling_regions.hg38.interval_list') )
        .set{ target_regions }

    Channel
        .empty()
        .set{ target_regions_bed }
} else if( params.seq_protocol == "WES" ) {
    Channel
        .value( file('references/hg38/wxs_exons_gencode_v39_autosome_sex_chroms.hg38.bed') )
        .set{ target_regions_bed }

    Channel
        .value( file('references/hg38/wxs_exons_gencode_v39_autosome_sex_chroms.hg38.interval_list') )
        .set{ target_regions }
} else {
    exit 1, "This run command cannot be executed. The '--seq_protocol' must be set to either 'WGS' for whole-genome or 'WES' for whole-exome."
}

Channel
	.fromPath( 'references/hg38/wgs_calling_regions.hg38.interval_list' )
	.set{ gatk_bundle_wgs_interval_list }

Channel
	.fromPath( 'references/hg38/Homo_sapiens_assembly38_autosome.interval_list' )
	.set{ autosome_chromosome_list }

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

log.info ''
log.info '################################################'
log.info ''
log.info "           .------------------------.           "
log.info "          |    .-..-. .--. .---.     |          "
log.info "          |    : `' :: .--': .; :    |          "
log.info "          |    : .. :: : _ :  _.'    |          "
log.info "          |    : :; :: :; :: :       |          "
log.info "          |    :_;:_;`.__.':_;       |          "
log.info "          |   ,-. .--.  .--.  .--.   |          "
log.info "          | .'  :: ,. :: ,. :: ,. :  |          "
log.info "          |   : :: :: :: :: :: :: :  |          "
log.info "          |   : :: :; :: :; :: :; :  |          "
log.info "          |   :_;`.__.'`.__.'`.__.'  |          "
log.info "           .________________________.           "
log.info ''
log.info "░█▀█░█▀▄░█▀▀░█▀█░█▀▄░█▀█░█▀▀░█▀▀░█▀▀░█▀▀░▀█▀░█▀█░█▀▀"
log.info "░█▀▀░█▀▄░█▀▀░█▀▀░█▀▄░█░█░█░░░█▀▀░▀▀█░▀▀█░░█░░█░█░█░█"
log.info "░▀░░░▀░▀░▀▀▀░▀░░░▀░▀░▀▀▀░▀▀▀░▀▀▀░▀▀▀░▀▀▀░▀▀▀░▀░▀░▀▀▀"
log.info ''
log.info "~~~ Launch Time ~~~"
log.info ''
log.info ' ${workflowTimestamp}'
log.info ''
log.info "~~~ Input Directory ~~~"
log.info ''
log.info ' ${params.input_dir}'
log.info ''
log.info "~~~ Output Directory ~~~"
log.info ''
log.info ' ${params.output_dir}'
log.info ''
log.info "~~~ Run Report File ~~~"
log.info ''
log.info ' nextflow_report.${params.run_id}.html'
log.info ''
log.info "~~~ Sequencing Protocol ~~~"
log.info ''
log.info " ${params.seq_protocol}"
log.info ''
log.info '################################################'
log.info ''

// if input files are BAMs, set the up channels for them to go through the pipeline or straight to BAM QC process
if( params.input_format == "bam" ) {
	Channel
		.fromPath( "${params.input_dir}/*.bam" )
		.ifEmpty{ error "BAM format specified but cannot find files with .bam extension in input directory" }
		.into{ input_mapped_bams; 
		       input_mapped_bams_forQaulimap }
} else {
	Channel
		.empty()
		.into{ input_mapped_bams;
			   input_mapped_bams_forQaulimap }
}

// If input files are FASTQs, set channel up for both R1 and R2 reads then merge into single channel
if( params.input_format == "fastq" ) {
	Channel
		.fromPath( "${params.input_dir}/*R{1,2}*.f*q*")
		.collect()
		.ifEmpty{ error "FASTQ format specified but cannot find files with expected R1/R2 naming convention, check test samples for example" }
		.set{ input_fastqs }
} else {
	Channel
		.empty()
		.set{ input_fastqs }
}

// Depending on if the input FASTQs needed be lane merged before being gathered
if( params.input_format == "fastq" & params.lane_split == "yes" ) {
	input_fastqs_forMerging = input_fastqs
}
else {
	input_fastqs_forMerging = Channel.empty()
}

// Lane-Split FASTQ Merge ~ for all input lane-split FASTQs, merge into single R1/R2 FASTQ file without altering input
process mergeLaneSplitFastqs_mergelane {
	publishDir "${params.output_dir}/preprocessing/", mode: 'symlink'

	input:
	path split_fastqs from input_fastqs_forMerging

	output:
	path lane_merged_input_fastqs into lane_merged_fastq_dir

	when:
	params.input_format == "fastq" & params.lane_split == "yes"

	script:
	lane_merged_input_fastqs = "laneMergedFastqs"
	"""
	lane_split_merger.sh \
	"${params.input_dir}" \
	. \
	"${lane_merged_input_fastqs}"
	"""
}

// If input FASTQs were lane merged, set as input for FASTQ gathering
if( params.input_format == "fastq" & params.lane_split == "yes" ) {
	fastqs_forGathering = lane_merged_fastq_dir
}
else {
	fastqs_forGathering = input_fastqs
}

// FASTQ Pair Gatherer ~ properly pair all input FASTQs and create sample sheet
process gatherInputFastqs_fastqgatherer {
	publishDir "${params.output_dir}/preprocessing/", mode: 'copy', pattern: '*.{txt}'

	input:
	path fastqs_forGathering

	output:
	path run_fastq_samplesheet into input_fastq_sample_sheet

	when:
	params.input_format == "fastq"

	script:
	run_fastq_samplesheet = "${params.run_id}.fastq.samplesheet.txt"
	"""
	fastq_pair_gatherer.pl \
	"${fastqs_forGathering}" \
	"${run_fastq_samplesheet}"
	"""
}

// If input files are FASTQs, read the input FASTQ sample sheet to set correct FASTQ pairs,
// then set channel up for both R1 and R2 reads then merge into single channel
if( params.input_format == "fastq" & params.lane_split == "yes" ) {
	input_fastq_sample_sheet.splitCsv( header: true, sep: '\t' )
						    .map{ row -> sample_id = "${row.sample_id}"
						                 input_R1_fastq = "${row.read_1}"
						                 input_R2_fastq = "${row.read_2}"
						          return[ "${sample_id}",
						                  file("${params.output_dir}/preprocessing/laneMergedFastqs/${input_R1_fastq}"),
						                  file("${params.output_dir}/preprocessing/laneMergedFastqs/${input_R2_fastq}") ] }
						    .set{ paired_input_fastqs }
} else if( params.input_format == "fastq" & params.lane_split == "no" & params.skip_trimming == "yes" ) {
	input_fastq_sample_sheet.splitCsv( header: true, sep: '\t' )
						    .map{ row -> sample_id = "${row.sample_id}"
						                 input_R1_fastq = "${row.read_1}"
						                 input_R2_fastq = "${row.read_2}"
						          return[ "${sample_id}",
						                  file("${params.input_dir}/${input_R1_fastq}"),
						                  file("${params.input_dir}/${input_R2_fastq}") ] }
						    .into{ paired_input_fastqs_forFastqc;
						           paired_input_fastqs_forAlignment }
} else if( params.input_format == "fastq" & params.lane_split == "no" & params.skip_trimming == "no" ) {
	input_fastq_sample_sheet.splitCsv( header: true, sep: '\t' )
						    .map{ row -> sample_id = "${row.sample_id}"
						                 input_R1_fastq = "${row.read_1}"
						                 input_R2_fastq = "${row.read_2}"
						          return[ "${sample_id}",
						                  file("${params.input_dir}/${input_R1_fastq}"),
						                  file("${params.input_dir}/${input_R2_fastq}") ] }
						    .set{ paired_input_fastqs }
} else {
	Channel
		.empty()
		.set{ paired_input_fastqs }
}

// GATK RevertSam ~ convert input mapped BAM files to unmapped BAM files
process revertMappedBam_gatk {
	tag "${sample_id}"

	input:
	path bam_mapped from input_mapped_bams

	output:
	path bam_unmapped into unmapped_bams

	when:
	params.input_format == "bam"

	script:
	bam_unmapped = "${bam_mapped}".replaceFirst(/\..*bam/, ".unmapped.bam")
	sample_id = "${bam_unmapped}".replaceFirst(/\.unmapped\.bam/, "")
	"""
	gatk RevertSam \
	--java-options "-Xmx${task.memory.toGiga() - 2}G -Djava.io.tmpdir=." \
	--VERBOSITY ERROR \
	--VALIDATION_STRINGENCY LENIENT \
	--MAX_RECORDS_IN_RAM 4000000 \
	--TMP_DIR . \
	--SANITIZE true \
	--ATTRIBUTE_TO_CLEAR XT \
	--ATTRIBUTE_TO_CLEAR XN \
	--ATTRIBUTE_TO_CLEAR OC \
	--ATTRIBUTE_TO_CLEAR OP \
	--INPUT "${bam_mapped}" \
	--OUTPUT "${bam_unmapped}"
	"""
}

// biobambam bamtofastq ~ convert unmapped BAM files to paired FASTQ files
process bamToFastq_biobambam {
	tag "${sample_id}"

	input:
	path bam_unmapped from unmapped_bams

	output:
	tuple val(sample_id), path(fastq_R1), path(fastq_R2) into converted_fastqs_forTrimming

	when:
	params.input_format == "bam"

	script:
	sample_id = "${bam_unmapped}".replaceFirst(/\.unmapped\.bam/, "")
	fastq_R1 = "${bam_unmapped}".replaceFirst(/\.unmapped\.bam/, "_R1.fastq.gz")
	fastq_R2 = "${bam_unmapped}".replaceFirst(/\.unmapped\.bam/, "_R2.fastq.gz")
	"""
	bamtofastq \
	filename="${bam_unmapped}" \
	F="${fastq_R1}" \
	F2="${fastq_R2}" \
	gz=1
	"""
}

// Depending on which input data type was used, set an input variable for the Trimmomatic process
if( params.input_format == "bam" & params.skip_trimming == "no" ) {
	input_fastqs_forTrimming = converted_fastqs_forTrimming
} else if( params.input_format == "fastq" & params.skip_trimming == "no" ) {
	input_fastqs_forTrimming = paired_input_fastqs
} else if( params.input_format == "fastq" & params.skip_trimming == "yes" ) {
	input_fastqs_forTrimming = Channel.empty()
}

// Trimmomatic ~ trim low quality bases and clip adapters from reads
process fastqTrimming_trimmomatic {
	publishDir "${params.output_dir}/preprocessing/trimLogs", mode: 'copy', pattern: '*.{log}'
	tag "${sample_id}"

	input:
	tuple val(sample_id), path(input_R1_fastqs), path(input_R2_fastqs), path(trimmomatic_contaminants) from input_fastqs_forTrimming.combine(trimmomatic_contaminants)

	output:
	tuple val(sample_id), path(fastq_R1_trimmed), path(fastq_R2_trimmed) into trimmed_fastqs_forFastqc, trimmed_fastqs_forAlignment
	path fastq_trim_log

	when:
	params.skip_trimming == "no"

	script:
	fastq_R1_trimmed = "${sample_id}_R1_trim.fastq.gz"
	fastq_R2_trimmed = "${sample_id}_R2_trim.fastq.gz"
	fastq_R1_unpaired = "${sample_id}_R1_unpaired.fastq.gz"
	fastq_R2_unpaired = "${sample_id}_R2_unpaired.fastq.gz"
	fastq_trim_log = "${sample_id}.trim.log"
	"""
	trimmomatic PE \
	-threads ${task.cpus} \
	"${input_R1_fastqs}" \
	"${input_R2_fastqs}" \
	"${fastq_R1_trimmed}" \
	"${fastq_R1_unpaired}" \
	"${fastq_R2_trimmed}" \
	"${fastq_R2_unpaired}" \
	ILLUMINACLIP:${trimmomatic_contaminants}:2:30:10:1:true \
	TRAILING:5 \
	SLIDINGWINDOW:4:15 \
	MINLEN:${params.trimmomatic_min_len} \
	2> "${fastq_trim_log}"
	"""
}

// Depending on if the FASTQs were trimmed or not, set the input for FastQC and alignment
if( params.skip_trimming == "no" ) {
	fastqs_forFastqc = trimmed_fastqs_forFastqc
	fastqs_forAlignment = trimmed_fastqs_forAlignment
} else if( params.skip_trimming == "yes" ) {
	fastqs_forFastqc = paired_input_fastqs_forFastqc
	fastqs_forAlignment = paired_input_fastqs_forAlignment
}

// FastQC ~ generate sequence quality metrics for input FASTQ files
process fastqQualityControlMetrics_fastqc {
	publishDir "${params.output_dir}/preprocessing/fastqc", mode: 'copy'
	tag "${sample_id}"

	input:
	tuple val(sample_id), path(fastq_R1), path(fastq_R2) from fastqs_forFastqc

	output:
	tuple path(fastqc_R1_html), path(fastqc_R2_html)
	tuple path(fastqc_R1_zip), path(fastqc_R2_zip)

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

// BWA MEM / Sambamba ~ align trimmed FASTQ files to reference genome to produce BAM file
process alignment_bwa {
	tag "${sample_id}"

	input:
	tuple val(sample_id), path(fastq_R1), path(fastq_R2), path(bwa_reference_dir) from fastqs_forAlignment.combine(bwa_reference_dir)

	output:
	path bam_aligned into aligned_bams
	tuple val(sample_id), path(bam_aligned) into aligned_bam_forFlagstats

	script:
	bam_aligned = "${sample_id}.bam"
	"""
	bwa mem \
	-M \
	-K 100000000 \
	-v 1 \
	-t ${task.cpus - 1} \
	-R '@RG\\tID:${sample_id}\\tSM:${sample_id}\\tLB:${sample_id}\\tPL:ILLUMINA' \
	"${bwa_reference_dir}/Homo_sapiens_assembly38.fasta" \
	"${fastq_R1}" "${fastq_R2}" \
	| \
	sambamba view \
	--sam-input \
	--nthreads=${task.cpus - 1} \
	--filter='mapping_quality>=10' \
	--format=bam \
	--compression-level=0 \
	/dev/stdin \
	| \
	sambamba sort \
	--nthreads=${task.cpus - 1} \
	--tmpdir=. \
	--memory-limit=24GB \
	--sort-by-name \
	--out=${bam_aligned} \
	/dev/stdin
	"""
}

// Sambamba flagstat ~ generate read metrics after alignment
process postAlignmentFlagstats_sambamba {
	publishDir "${params.output_dir}/preprocessing/alignmentFlagstats", mode: 'copy', pattern: "*${bam_flagstat_log}"
	tag "${sample_id}"

	input:
	tuple val(sample_id), path(bam_aligned) from aligned_bam_forFlagstats

	output:
	path bam_flagstat_log

	script:
	bam_flagstat_log = "${sample_id}.alignment.flagstat.log"
	"""
	sambamba flagstat \
	--nthreads=${task.cpus} \
	"${bam_aligned}" > "${bam_flagstat_log}"
	"""
}

// GATK FixMateInformation / SortSam ~ veryify/fix mate-pair information and sort output BAM by coordinate
process fixMateInformationAndSort_gatk {
	tag "${bam_aligned.baseName}"

	input:
	path bam_aligned from aligned_bams

	output:
	path bam_fixed_mate into fixed_mate_bams

	script:
	bam_fixed_mate_unsorted = "${bam_aligned}".replaceFirst(/\.bam/, ".unsorted.fixedmate.bam")
	bam_fixed_mate = "${bam_aligned}".replaceFirst(/\.bam/, ".fixedmate.bam")
	"""
	gatk FixMateInformation \
	--java-options "-Xmx16G -Djava.io.tmpdir=. -XX:ParallelGCThreads=1" \
	--VERBOSITY ERROR \
	--VALIDATION_STRINGENCY SILENT \
	--ADD_MATE_CIGAR true \
	--MAX_RECORDS_IN_RAM 2000000 \
	--ASSUME_SORTED true \
	--TMP_DIR . \
	--INPUT "${bam_aligned}" \
	--OUTPUT "${bam_fixed_mate_unsorted}"

	gatk SortSam \
	--java-options "-Xmx16G -Djava.io.tmpdir=. -XX:ParallelGCThreads=1" \
	--VERBOSITY ERROR \
	--TMP_DIR . \
	--SORT_ORDER coordinate \
	--INPUT "${bam_fixed_mate_unsorted}" \
	--OUTPUT "${bam_fixed_mate}"
	"""
}

// Sambamba markdup ~ mark duplicate alignments, remove them, and create BAM index
process markDuplicatesAndIndex_sambamba {
	publishDir "${params.output_dir}/preprocessing/markdupFlagstats", mode: 'copy', pattern: '*.{log}'
	tag "${sample_id}"

	input:
	path bam_fixed_mate from fixed_mate_bams

	output:
	tuple val(sample_id), path(bam_marked_dup) into marked_dup_bams_forRealignment, marked_dup_bams_forDownsampleBam, marked_dup_bams_forApplyBqsr
	path bam_marked_dup_index
	path markdup_output_log
	path bam_markdup_flagstat_log

	script:
	sample_id = "${bam_fixed_mate}".replaceFirst(/\.fixedmate\.bam/, "")
	bam_marked_dup = "${sample_id}.markdup.bam"
	bam_marked_dup_index = "${bam_marked_dup}.bai"
	markdup_output_log = "${sample_id}.markdup.log"
	bam_markdup_flagstat_log = "${sample_id}.markdup.flagstat.log"
	"""
	sambamba markdup \
	--remove-duplicates \
	--nthreads ${task.cpus} \
	--hash-table-size 1000000 \
	--overflow-list-size 1000000 \
	--tmpdir . \
	"${bam_fixed_mate}" \
	"${bam_marked_dup}" \
	2> "${markdup_output_log}"

	sambamba flagstat \
	--nthreads ${task.cpus} \
	"${bam_marked_dup}" > "${bam_markdup_flagstat_log}"

	sambamba index \
	--nthreads ${task.cpus} \
	"${bam_marked_dup}" "${bam_marked_dup_index}"
	"""	
}

reference_genome_fasta_forRealignment.combine( reference_genome_fasta_index_forRealignment )
    .combine( reference_genome_fasta_dict_forRealignment )
    .set{ reference_genome_and_targets_bundle_forRealignment }

// ABRA2 ~ local and global realignment for improvement of InDel calling in exome data
process localAndGlobalRealignment_abra2 {
    publishDir "${params.output_dir}/preprocessing/realignment", mode: 'copy', pattern: '*.{log}'
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(bam_marked_dup), path(reference_genome_fasta_forRealignment), path(reference_genome_fasta_index_forRealignment), path(reference_genome_fasta_dict_forRealignment) from marked_dup_bams_forRealignment.combine(reference_genome_and_targets_bundle_forRealignment)
    path target_bed from target_regions_bed

    output:
    tuple val(sample_id), path(bam_marked_dup_realigned) into marked_dup_realigned_bams_forDownsampleBam, marked_dup_realigned_bams_forApplyBqsr
    path abra_log

    when:
    params.seq_protocol == "WES"

    script:
    bam_marked_dup_realigned = "${bam_marked_dup}".replaceFirst(/\.bam/, ".realign.bam")
    abra_log = "${sample_id}.abra.realign.log"
    """
    java -jar -Xmx16G \$ABRA2JAR \
    --in "${bam_marked_dup}" \
    --out "${bam_marked_dup_realigned}" \
    --ref "${reference_genome_fasta_forRealignment}" \
    --targets "${target_bed}" \
    --threads ${task.cpus} \
    --tmpdir . \
    > "${abra_log}"
    """
}

if( params.seq_protocol == "WGS" ) {
    bams_forDownsampleBam = marked_dup_bams_forDownsampleBam
} else if( params.seq_protocol == "WES" ) {
    bams_forDownsampleBam = marked_dup_realigned_bams_forDownsampleBam
}

// GATK DownsampleSam ~ downsample BAM file to use random subset for generating BSQR table
process downsampleBam_gatk {
	tag "${sample_id}"

	input:
	tuple val(sample_id), path(bam_for_downsample) from bams_forDownsampleBam

	output:
	tuple val(sample_id), path(bam_downsampled) into downsampled_bams

	script:
	bam_downsampled = "${bam_for_downsample}".replaceFirst(/\.bam/, ".downsampled.bam")
	"""
	gatk DownsampleSam \
	--java-options "-Xmx${task.memory.toGiga() - 2}G -Djava.io.tmpdir=." \
	--VERBOSITY ERROR \
	--MAX_RECORDS_IN_RAM 4000000 \
	--TMP_DIR . \
	--STRATEGY ConstantMemory \
	--RANDOM_SEED 1000 \
	--CREATE_INDEX \
	--VALIDATION_STRINGENCY SILENT \
	--PROBABILITY 0.1 \
	--INPUT "${bam_for_downsample}" \
	--OUTPUT "${bam_downsampled}"
	"""
}

// Combine all needed GATK bundle files and reference FASTA into one channel for use in GATK BaseRecalibrator process
gatk_bundle_mills_1000G.combine( gatk_bundle_mills_1000G_index )
	.combine( gatk_bundle_known_indels )
	.combine( gatk_bundle_known_indels_index )
	.combine( gatk_bundle_dbsnp138 )
	.combine( gatk_bundle_dbsnp138_index )
	.set{ gatk_reference_bundle }

reference_genome_fasta_forBaseRecalibrator.combine( reference_genome_fasta_index_forBaseRecalibrator )
	.combine( reference_genome_fasta_dict_forBaseRecalibrator )
	.set{ reference_genome_bundle_forBaseRecalibrator }

// Combine the the input BAM, GATK bundle, and reference FASTA files into one channel
downsampled_bams.combine( reference_genome_bundle_forBaseRecalibrator )
	.combine( gatk_reference_bundle )
	.set{ input_and_reference_files_forBaseRecalibrator }

// GATK BaseRecalibrator ~ generate base quality score recalibration table based on covariates
process baseRecalibrator_gatk {
	tag "${sample_id}"

	input:
	tuple val(sample_id), path(bam_downsampled), path(reference_genome_fasta_forBaseRecalibrator), path(reference_genome_fasta_index_forBaseRecalibrator), path(reference_genome_fasta_dict_forBaseRecalibrator), path(gatk_bundle_mills_1000G), path(gatk_bundle_mills_1000G_index), path(gatk_bundle_known_indels), path(gatk_bundle_known_indels_index), path(gatk_bundle_dbsnp138), path(gatk_bundle_dbsnp138_index) from input_and_reference_files_forBaseRecalibrator
	path targets_list from target_regions

	output:
	tuple val(sample_id), path(bqsr_table) into base_quality_score_recalibration_data

	script:
	bqsr_table = "${sample_id}.recaldata.table"
	"""
	gatk BaseRecalibrator \
	--java-options "-Xmx${task.memory.toGiga() - 2}G -Djava.io.tmpdir=." \
	--verbosity ERROR \
	--tmp-dir . \
	--read-filter GoodCigarReadFilter \
	--reference "${reference_genome_fasta_forBaseRecalibrator}" \
	--intervals "${targets_list}" \
	--input "${bam_downsampled}" \
	--output "${bqsr_table}" \
	--known-sites "${gatk_bundle_mills_1000G}" \
	--known-sites "${gatk_bundle_known_indels}" \
	--known-sites "${gatk_bundle_dbsnp138}"
	"""
}

if( params.seq_protocol == "WGS" ) {
    bams_forApplyBqsr = marked_dup_bams_forApplyBqsr
} else if( params.seq_protocol == "WES" ) {
    bams_forApplyBqsr = marked_dup_realigned_bams_forApplyBqsr
}

// Create additional channel for the reference FASTA to be used in GATK ApplyBQSR process
reference_genome_fasta_forApplyBqsr.combine( reference_genome_fasta_index_forApplyBqsr )
	.combine( reference_genome_fasta_dict_forApplyBqsr )
	.set{ reference_genome_bundle_forApplyBqsr }

// First join the input BAM files with their respective BQSR recalibration table, then combine that with the
// reference FASTA files into one channel
bams_forApplyBqsr.join( base_quality_score_recalibration_data )
	.set{ bams_and_bqsr_tables }

bams_and_bqsr_tables.combine( reference_genome_bundle_forApplyBqsr )
	.set{ input_and_reference_files_forApplyBqsr }

// GATK ApplyBQSR ~ apply base quality score recalibration using generated table
process applyBqsr_gatk {
	publishDir "${params.output_dir}/preprocessing/finalPreprocessedBams", mode: 'copy', pattern: '*.{final.bam,bai}'
	tag "${sample_id}"

	input:
	tuple val(sample_id), path(bam_for_bqsr), path(bqsr_table), path(reference_genome_fasta_forApplyBqsr), path(reference_genome_fasta_index_forApplyBqsr), path(reference_genome_fasta_dict_forApplyBqsr) from input_and_reference_files_forApplyBqsr

	output:
	tuple val(sample_id), path(bam_preprocessed_final) into final_preprocessed_bams_forCollectWgsMetrics, final_preprocessed_bams_forCollectGcBiasMetrics, final_preprocessed_bams_forCollectHsMetrics
	path bam_preprocessed_final_index

	script:
	bam_preprocessed_final = "${bam_for_bqsr}".replaceFirst(/\..*\.bam/, ".final.bam")
	bam_preprocessed_final_index = "${bam_preprocessed_final}".replaceFirst(/\.bam$/, ".bai")
	"""
	gatk ApplyBQSR \
	--java-options "-Xmx${task.memory.toGiga() - 2}G -Djava.io.tmpdir=." \
	--verbosity ERROR \
	--tmp-dir . \
	--read-filter GoodCigarReadFilter \
	--reference "${reference_genome_fasta_forApplyBqsr}" \
	--input "${bam_for_bqsr}" \
	--output "${bam_preprocessed_final}" \
	--bqsr-recal-file "${bqsr_table}"
	"""
}

// Create additional channel for the reference FASTA and autosome chromosome only interval list to be used in GATK CollectWgsMetrics process
reference_genome_fasta_forCollectWgsMetrics.combine( reference_genome_fasta_index_forCollectWgsMetrics )
	.combine( reference_genome_fasta_dict_forCollectWgsMetrics )
	.combine( autosome_chromosome_list )
	.set{ reference_genome_bundle_forCollectWgsMetrics }

// GATK CollectWgsMetrics ~ generate covearge and performance metrics from final BAM
process collectWgsMetrics_gatk {
	publishDir "${params.output_dir}/preprocessing/coverageMetrics", mode: 'copy'
	tag "${sample_id}"

	input:
	tuple val(sample_id), path(bam_preprocessed_final), path(reference_genome_fasta_forCollectWgsMetrics), path(reference_genome_fasta_index_forCollectWgsMetrics), path(reference_genome_fasta_dict_forCollectWgsMetrics), path(autosome_chromosome_list) from final_preprocessed_bams_forCollectWgsMetrics.combine( reference_genome_bundle_forCollectWgsMetrics)

	output:
	path coverage_metrics

	when:
	params.seq_protocol == "WGS"

	script:
	coverage_metrics = "${sample_id}.coverage.metrics.txt"
	"""
	gatk CollectWgsMetrics \
	--java-options "-Xmx${task.memory.toGiga() - 2}G -Djava.io.tmpdir=." \
	--VERBOSITY ERROR \
	--TMP_DIR . \
	--INCLUDE_BQ_HISTOGRAM \
	--MINIMUM_BASE_QUALITY 20 \
	--MINIMUM_MAPPING_QUALITY 20 \
	--REFERENCE_SEQUENCE "${reference_genome_fasta_forCollectWgsMetrics}" \
	--INTERVALS "${autosome_chromosome_list}" \
	--INPUT "${bam_preprocessed_final}" \
	--OUTPUT "${coverage_metrics}"
	"""
}

// Create additional channel for the reference FASTA and interfal list to be used in GATK CollectWgsMetrics process
reference_genome_fasta_forCollectGcBiasMetrics.combine( reference_genome_fasta_index_forCollectGcBiasMetrics )
	.combine( reference_genome_fasta_dict_forCollectGcBiasMetrics )
	.set{ reference_genome_bundle_forCollectGcBiasMetrics }

// GATK CollectGcBiasMetrics ~ generate GC content bias in reads in final BAM
process collectGcBiasMetrics_gatk {
	publishDir "${params.output_dir}/preprocessing/gcBiasMetrics", mode: 'copy'
	tag "${sample_id}"

	input:
	tuple val(sample_id), path(bam_preprocessed_final), path(reference_genome_fasta_forCollectGcBiasMetrics), path(reference_genome_fasta_index_forCollectGcBiasMetrics), path(reference_genome_fasta_dict_forCollectGcBiasMetrics) from final_preprocessed_bams_forCollectGcBiasMetrics.combine(reference_genome_bundle_forCollectGcBiasMetrics)

	output:
	path gc_bias_metrics
	path gc_bias_chart
	path gc_bias_summary

	script:
	gc_bias_metrics = "${sample_id}.gcbias.metrics.txt"
	gc_bias_chart = "${sample_id}.gcbias.metrics.pdf"
	gc_bias_summary = "${sample_id}.gcbias.summary.txt"
	"""
	gatk CollectGcBiasMetrics \
	--java-options "-Xmx${task.memory.toGiga() - 2}G -Djava.io.tmpdir=." \
	--VERBOSITY ERROR \
	--TMP_DIR . \
	--REFERENCE_SEQUENCE "${reference_genome_fasta_forCollectGcBiasMetrics}" \
	--INPUT "${bam_preprocessed_final}" \
	--OUTPUT "${gc_bias_metrics}" \
	--CHART_OUTPUT "${gc_bias_chart}" \
	--SUMMARY_OUTPUT "${gc_bias_summary}"
	"""
}

// Create additional channel for the reference FASTA to be used in GATK CollectHsMetrics process
reference_genome_fasta_forCollectHsMetrics.combine( reference_genome_fasta_index_forCollectHsMetrics )
	.combine( reference_genome_fasta_dict_forCollectHsMetrics )
	.set{ reference_genome_bundle_forCollectHsMetrics }

// GATK CollectHsMetrics ~ generate exome capture coverage data in final BAM
process collectHsMetrics_gatk {
    publishDir "${params.output_dir}/preprocessing/hsMetrics", mode: 'copy'
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(bam_preprocessed_final), path(reference_genome_fasta_forCollectHsMetrics), path(reference_genome_fasta_index_forCollectHsMetrics), path(reference_genome_fasta_dict_forCollectHsMetrics) from final_preprocessed_bams_forCollectHsMetrics.combine(reference_genome_bundle_forCollectHsMetrics)
    path targets_list from target_regions

    output:
    path hs_metrics

    when:
    params.seq_protocol == "WES"

    script:
    hs_metrics = "${sample_id}.hs.metrics.txt"
    """
    gatk CollectHsMetrics \
    --java-options "-Xmx${task.memory.toGiga() - 2}G -Djava.io.tmpdir=." \
    --VERBOSITY ERROR \
    --TMP_DIR . \
    --COVERAGE_CAP 500 \
    --REFERENCE_SEQUENCE "${reference_genome_fasta_forCollectHsMetrics}" \
    --BAIT_INTERVALS "${targets_list}" \
    --TARGET_INTERVALS "${targets_list}" \
    --INPUT "${bam_preprocessed_final}" \
    --OUTPUT "${hs_metrics}" \
    """
}
