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




params.output_dir = "output"

// ################################################ \\
//                                                  \\
// ~~~~~~~~~~~~~~~~ PRE-PROCESSING ~~~~~~~~~~~~~~~~ \\
//                                                  \\
// ################################################ \\

input_bams = Channel.fromPath( 'testData/*.fastq.gz' )

// FastQC ~ generate sequence quality metrics for input BAM files
process fastqc {
	publishDir "${params.output_dir}/preprocessing/fastqc", mode: 'copy'

	input:
	path bam from input_bams

	output:
	path fastqc_html into fastqc_reports
	path fastqc_zip

	script:
	fastqc_html = "${bam}".replaceFirst(/.fastq.gz$/, "_fastqc.html")
	fastqc_zip = "${bam}".replaceFirst(/.fastq.gz$/, "_fastqc.zip")
	"""
	fastqc -o . "${bam}"
	"""
}








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