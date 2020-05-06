//  |\        /|   /---------\  /---------\
//  | \      / |  |             |         |
//  |  \    /  |  |   /------\  |---------/
//  |   \  /   |  |          |  |
//  |    \/    |   \---------/  |

input_fastqs = Channel.fromPath('~/morganLab/mgp1000/development/testData/*.fastqz.gz')
params.output_dir = "output"

process fastqc {
	publishDir "${params.output_dir}"
	echo true

	input:
	file(fastq) from input_fastqs

	output:
	file(output_html)
	file(output_zip)

	script:
	"""
	fastqc -o . "${fastq}"
	"""
}