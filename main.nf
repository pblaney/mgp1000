//  ################################################################################
//  |\       /|   /---------\  /---------\    /|   /------\    /------\    /------\
//  | \     / |  |             |         |   / |  |        |  |        |  |        |
//  |  \   /  |  |     |----\  |---------/     |  |        |  |        |  |        |
//  |   \ /   |  |          |  |               |  |        |  |        |  |        |
//  |         |   \---------/  |               |   \------/    \------/    \------/
//  ################################################################################

params.output_dir = "output"

input_fastqs = Channel.fromPath( 'testData/*.fastq.gz' )
process fastqc {
	publishDir "${params.output_dir}/fastqc", mode: 'copy'
	containerOptions '-v $PWD/testData:/home/data'

	input:
	file fastq from input_fastqs

	output:
	file "*_fastqc.{zip,html}" into fastqc_reports

	script:
	"""
	fastqc "/home/data/${fastq}"
	"""
}

process multiqc {
	publishDir "${params.output_dir}/multiqc", mode: 'copy'

	input:
	file fastqc_dir from fastqc_reports.collect()

	output:
	file "multiqc_report.html"
	file "multiqc_data"

	script:
	"""
	multiqc "${fastqc_dir}"
	"""
}