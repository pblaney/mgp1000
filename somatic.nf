// Myeloma Genome Project 1000
// Comprehensive pipeline for analysis of matched T/N Multiple Myeloma WGS data
// https://github.com/pblaney/mgp1000

// This portion of the pipeline is used for somatic variant analysis of matched tumor/normal WGS samples.
// It is designed to be run with BAMs that were genereated via the Preprocessing step of this pipeline.

import java.text.SimpleDateFormat;
def workflowTimestamp = "${workflow.start.format('MM-dd-yyyy HH:mm')}"

log.info ''
log.info '##### Myeloma Genome Project 1000 Pipeline #####'
log.info '################################################'
log.info '~~~~~~~~~~~ SOMATIC VARIANT ANALYSIS ~~~~~~~~~~~'
log.info '################################################'
log.info ''
log.info "~~~ Launch Time ~~~		${workflowTimestamp}"
log.info ''
log.info "~~~ Output Directory ~~~ 	${workflow.projectDir}/output/somatic"
log.info ''
log.info "~~~ Run Report File ~~~ 	nextflow_report.somatic_${params.run_id}.html"
log.info ''
log.info '################################################'
log.info ''

def helpMessage() {
	log.info"""

	Usage Example:

		nextflow run somatic.nf -bg -resume --run_id batch1 --sample_sheet samplesheet.csv --singularity_module singularity/3.1 --email someperson@gmail.com -profile somatic 

	Mandatory Arguments:
    	--run_id                       [str]  Unique identifier for pipeline run
    	--sample_sheet                 [str]  CSV file containing the list of samples where the first column designates the file name of the
    	                                      normal sample, the second column for the file name of the matched tumor sample, example of the
    	                                      format for this file is in the testSamples directory
		-profile                       [str]  Configuration profile to use, each profile described in nextflow.config file
		                                      Currently available: preprocessing, germline, somatic

	Main Options:
		-bg                           [flag]  Runs the pipeline processes in the background, this option should be included if deploying
		                                      pipeline with real data set so processes will not be cut if user disconnects from deployment
		                                      environment
		-resume                       [flag]  Successfully completed tasks are cached so that if the pipeline stops prematurely the
		                                      previously completed tasks are skipped while maintaining their output
		--email                        [str]  Email address to send workflow completion/stoppage notification
		--singularity_module           [str]  Indicates the name of the Singularity software module to be loaded for use in the pipeline,
		                                      this option is not needed if Singularity is natively installed on the deployment environment
		--help                        [flag]  Prints this message

	################################################

	""".stripIndent()
}

// #################################################### \\
// ~~~~~~~~~~~~~ PARAMETER CONFIGURATION ~~~~~~~~~~~~~~ \\

// Declare the defaults for all pipeline parameters
params.output_dir = "output"
params.run_id = null
params.sample_sheet = null
params.help = null

// Print help message if requested
if( params.help ) exit 0, helpMessage()

// Print error messages if required parameters are not set
if( params.run_id == null ) exit 1, "The run command issued does not have the '--run_id' parameter set. Please set the '--run_id' parameter to a unique identifier for the run."

if( params.sample_sheet == null ) exit 1, "The run command issued does not have the '--sample_sheet' parameter set. Please set the '--sample_sheet' parameter to the path of the normal/tumor pair sample sheet CSV."

// Set channels for reference files
Channel
	.fromPath( 'references/hg38/Homo_sapiens_assembly38.fasta' )


Channel
	.fromPath( 'references/hg38/Homo_sapiens_assembly38.fasta.fai' )


Channel
	.fromPath( 'references/hg38/Homo_sapiens_assembly38.dict' )


// #################################################### \\
// ~~~~~~~~~~~~~~~~ PIPELINE PROCESSES ~~~~~~~~~~~~~~~~ \\

// Read user provided sample sheet to find Tumor/Normal paired sample BAM files
Channel
	.fromPath( params.sample_sheet )
	.ifEmpty{ error "No sample sheet provided, an example is given in the testSamples directory" }
	.splitCsv( header:true )
	.map{ row -> file("input/preprocessedBams/${row.normal}") }
	.into{ input_preprocessed_bams_for }