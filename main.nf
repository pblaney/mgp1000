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
trimmomatic_contaminants = Channel.fromPath('references/trimmomatic.fa')

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
	fastq_R1 = "${bam}".replaceFirst(/.bam$/, "_R1.fastq")
	fastq_R2 = "${bam}".replaceFirst(/.bam$/, "_R2.fastq")
	"""
	bamtofastq filename="${bam}" F="${fastq_R1}" F2="${fastq_R2}" gz=1
	"""
}

// Set the path to all input FASTQ files. If started with BAM files, use the converted FASTQs
if( params.input_format == "bam" ) {
	input_fastqs = converted_fastqs
}
else {
	input_R1_fastqs = Channel.fromPath( 'input/*R1.f*q*' )
	input_R2_fastqs = Channel.fromPath( 'input/*R2.f*q*' )
	input_fastqs = input_R1_fastqs.merge( input_R2_fastqs )
}

// Trimmomatic ~ trim low quality bases and clip adapters from reads
process trimmomatic {
	publishDir "${params.output_dir}/preprocessing/trimmomatic/trimmedFastqs", mode: 'copy', pattern: "*${fastq_R1_trimmed}"
	publishDir "${params.output_dir}/preprocessing/trimmomatic/trimmedFastqs", mode: 'copy', pattern: "*${fastq_R2_trimmed}"

	input:
	tuple path(input_R1_fastqs), path(input_R2_fastqs) from input_fastqs
	path trimmomatic_contaminants

	output:
	tuple path(fastq_R1_trimmed), path(fastq_R2_trimmed) into trimmed_fastqs

	script:
	fastq_R1_trimmed = "${input_R1_fastqs}".replaceFirst(/.fastq.gz$/, ".trim.fastq.gz")
	fastq_R2_trimmed = "${input_R2_fastqs}".replaceFirst(/.fastq.gz$/, ".trim.fastq.gz")
	fastq_R1_unpaired = "${input_R1_fastqs}".replaceFirst(/.fastq.gz$/, ".unpaired.fastq.gz")
	fastq_R2_unpaired = "${input_R2_fastqs}".replaceFirst(/.fastq.gz$/, ".unpaired.fastq.gz")
	"""
	trimmomatic PE -threads 8 \
	"${input_R1_fastqs}" "${input_R2_fastqs}" \
	"${fastq_R1_trimmed}" "${fastq_R1_unpaired}" \
	"${fastq_R2_trimmed}" "${fastq_R2_unpaired}" \
	ILLUMINACLIP:${trimmomatic_contaminants}:2:30:10:1:true TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:35
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







/*

// BWA MEM + Sambamba ~ align input BAM files to reference genome
process alignment {
	
}

*/



/* ###### WORK IN PROGRESS ######

process multiqc {
	publishDir "${params.output_dir}/multiqc", mode: 'copy'

	input:
	path reports from fastqc_reports.collect()

	output:
	path "multiqc_report.html"
	path "multiqc_data"

	script:
	"""
	multiqc "$PWD"
	"""
}

*/