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
	containerOptions '-v ~/morganLab/mgp1000/development/mgp1000/testData:/home/data'

	input:
	path fastq from input_fastqs

	output:
	path fastqc_html into fastqc_reports
	path fastqc_zip

	script:
	fastqc_html = "${fastq}".replaceFirst(/.fastq.gz$/, "_fastqc.html")
	fastqc_zip = "${fastq}".replaceFirst(/.fastq.gz$/, "_fastqc.zip")
	"""
	fastqc -o . "/home/data/${fastq}"
	"""
}

//process multiqc {
//	publishDir "${params.output_dir}/multiqc", mode: 'copy'
//
//	input:
//	path  from fastqc_reports.collect()
//
//	output:
//	path "multiqc_report.html"
//	path "multiqc_data"
//
//	script:
//	"""
//	multiqc .
//	"""
//}