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
params.input_format = "bam"



// ################################################ \\
//                                                  \\
// ~~~~~~~~~~~~~~~~ PRE-PROCESSING ~~~~~~~~~~~~~~~~ \\
//                                                  \\
// ################################################ \\

input_data = Channel.fromPath( 'input/' )


// biobambam ~ convert input BAM files to FASTQ files
process biobambam {
	publishDir "${params.out_dir}/preprocessing/biobambam", mode: 'copy'

	input:
	path bam from input_data

	output:
	tuple path(fastq_R1), path(fastq_R2) into input_fastqs 

	when:
	params.input_format != "bam"

	script:
	fastq_R1 = "${bam}".replaceFirst(/.bam$/, "_R1.fastq.gz")
	fastq_R2 = "${bam}".replaceFirst(/.bam$/, "_R2.fastq.gz")
	"""
	bamtofastq filename="${bam}" F="${fastq_R1}" F2="${fastq_R2}" gz=1
	"""
}




/*

input_fastqs = Channel.fromPath( 'testData/*.fastq.gz' )

// FastQC ~ generate sequence quality metrics for input FASTQ files
process fastqc {
	publishDir "${params.output_dir}/preprocessing/fastqc", mode: 'copy'

	input:
	path fastq from input_fastqs

	output:
	path fastqc_html into fastqc_reports
	path fastqc_zip

	script:
	fastqc_html = "${fastq}".replaceFirst(/.fastq.gz$/, "_fastqc.html")
	fastqc_zip = "${fastq}".replaceFirst(/.fastq.gz$/, "_fastqc.zip")
	"""
	fastqc --quite --noextract --outdir . "${fastq}"
	"""
}

// Trimmomatic ~ trim low quality bases and clip adapters from reads
process trimmomatic {
	
}


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