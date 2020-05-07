//  ################################################################################
//  |\       /|   /---------\  /---------\    /|   /------\    /------\    /------\
//  | \     / |  |             |         |   / |  |        |  |        |  |        |
//  |  \   /  |  |     |----\  |---------/     |  |        |  |        |  |        |
//  |   \ /   |  |          |  |               |  |        |  |        |  |        |
//  |         |   \---------/  |               |   \------/    \------/    \------/
//  ################################################################################

input_fastqs = Channel.fromPath( 'testData/*.fastq.gz' )
params.output_dir = "output"

process fastqc {
	tag "${fastq}"
	publishDir "${params.output_dir}/fastqc", mode: 'copy', overwrite: true
	echo true
	containerOptions '-v /Users/blanep01/morganLab/mgp1000/development/mgp1000/testData:/home/data'

	input:
	file fastq from input_fastqs

	output:
	file(output_html)
	file(output_zip)

	script:
	output_html = "${fastq}".replaceFirst(/.fastq.gz$/, "_fastqc.html")
	output_zip = "${fastq}".replaceFirst(/.fastq.gz$/, "_fastqc.zip")
	"""
	fastqc -o . "/home/data/${fastq}"
	"""
}