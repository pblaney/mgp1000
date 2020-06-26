// Myeloma Genome Project 1000
// Comprehensive pipeline for analysis of matched T/N Multiple Myeloma WGS data
// https://github.com/pblaney/mgp1000

// This portion of the pipeline is used for germline variant analysis of normal WGS samples.
// It can be run on any preprocessed BAM files, even if not run through the Preprocessing step of
// this pipeline.

log.info ''
log.info '##### Myeloma Genome Project 1000 Pipeline #####'
log.info '################################################'
log.info '~~~~~~~~~~~ GERMLINE VARIANT ANALYSIS ~~~~~~~~~~'
log.info '################################################'
log.info ''

// #################################################### \\
// ~~~~~~~~~~~~~ PARAMETER CONFIGURATION ~~~~~~~~~~~~~~ \\

// Declare the defaults for all pipeline parameters
params.output_dir = "output"

// #################################################### \\
// ~~~~~~~~~~~~~~~~ PIPELINE PROCESSES ~~~~~~~~~~~~~~~~ \\

// Set the channel for all input BAM files
Channel
	.fromPath( 'input/preprocessedBams/*.bam' )
	.into( 'input_preprocessed_bams_forTelomereLengthEstimation';
	       'input_preprocessed_bams_forHaplotypeCaller' )

// Telomerecat bam2length ~ estimate telomere length of sample
process telomereLengthEstimation_telomerecat {
	publishDir "${params.output_dir}/germline/telomereLength", mode: 'move'

	input:
	path bam_preprocessed from input_preprocessed_bams_forTelomereLengthEstimation

	output:
	path telomere_length_estimation

	script:
	telomere_length_estimation = "${bam_preprocessed}".replaceFirst(/\..*bam/, ".telomerelength.csv")
	"""
	telomerecat bam2length -p4 -v1 --output "${telomere_length_estimation}" "${bam_preprocessed}"
	"""
}