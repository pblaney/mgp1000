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

// Declare the defaults for all pipeline parameters
params.input_format = "bam"
params.output_dir = "output"

// Set path to reference files
trimmomatic_contaminants = Channel.fromPath( 'references/trimmomatic.fa' )
bwa_reference_genome = Channel.fromPath( 'references/hg38/bwa/genome.fa' )

// ################################################ \\
//                                                  \\
// ~~~~~~~~~~~~~~~~ PRE-PROCESSING ~~~~~~~~~~~~~~~~ \\
//                                                  \\
// ################################################ \\

// Set the path to all input BAM files
input_bams = Channel.fromPath( 'input/*.bam' )

// biobambam ~ convert any input BAM files to FASTQ files
process biobambam {
	publishDir "${params.output_dir}/preprocessing/biobambam", mode: 'copy'

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
process trimmomatic {
	publishDir "${params.output_dir}/preprocessing/trimmomatic/trimmedFastqs", mode: 'copy', pattern: "*${fastq_R1_trimmed}"
	publishDir "${params.output_dir}/preprocessing/trimmomatic/trimmedFastqs", mode: 'copy', pattern: "*${fastq_R2_trimmed}"
	publishDir "${params.output_dir}/preprocessing/trimmomatic/trimLogs", mode: 'copy', pattern: "*${fastq_trim_log}"

	input:
	tuple path(input_R1_fastqs), path(input_R2_fastqs) from input_fastqs
	path trimmomatic_contaminants

	output:
	tuple path(fastq_R1_trimmed), path(fastq_R2_trimmed) into trimmed_fastqs
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
process fastqc {
	publishDir "${params.output_dir}/preprocessing/fastqc", mode: 'copy'

	input:
	tuple path(fastq_R1), path(fastq_R2) from trimmed_fastqs

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
process alignment {
	publishDir "${params.output_dir}/preprocessing/alignment/rawBams", mode: 'copy', pattern: "*${bam_aligned}"
	publishDir "${params.output_dir}/preprocessing/alignment/flagstatLogs", mode: 'copy', pattern: "*${bam_flagstat_log}"

	input:
	tuple path(fastq_R1_trimmed), path(fastq_R2_trimmed) from trimmed_fastqs
	path bwa_reference_genome

	output:
	path bam_aligned into aligned_bams
	path bam_flagstat_log

	script:
	sample_id = "${fastq_R1_trimmed}".replaceFirst(/.trim.fastq.gz$/, "")
	bam_aligned = "${sample_id}.bam"
	bam_flagstat_log = "${sample_id}.alignFlagstat.log"
	"""
	bwa mem \
	-M -v 1 \
	-t 8 \
	-R '@RG\tID:${sample_id}\tSM:${sample_id}\tLB:${sample_id}\tPL:ILLUMINA' \
	"${bwa_reference_genome}" \
	"${fastq_R1_trimmed}" "${fastq_R2_trimmed}" \
	| \
	sambamba view \
	--sam-input \
	--nthreads=8 \
	--filter='mapping_quality>=10' \
	--format=bam \
	--compression-level=0 \
	/dev/stdin \
	| \
	sambamba sort \
	--nthreads=8 \
	--out="${bam_aligned}"

	sambamba flagstat "${bam_aligned}" > "${bam_flagstat_log}"
	"""
}

