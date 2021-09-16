// Myeloma Genome Project 1000
// Comprehensive pipeline for analysis of matched T/N Multiple Myeloma WGS data
// https://github.com/pblaney/mgp1000

// This portion of the pipeline is used for somatic variant analysis of matched tumor/normal WGS samples.
// It is designed to be run with BAMs that were genereated via the Preprocessing step of this pipeline.

import java.text.SimpleDateFormat;
def workflowTimestamp = "${workflow.start.format('MM-dd-yyyy HH:mm')}"

def helpMessage() {
	log.info"""

	Usage Example:

		nextflow run somatic.nf -bg -resume --run_id batch1 --sample_sheet samplesheet.csv --singularity_module singularity/3.1 --email someperson@gmail.com --mutect_ref_vcf_concatenated no --vep_ref_cached no -profile somatic 

	Mandatory Arguments:
    	--run_id                       [str]  Unique identifier for pipeline run
    	--sample_sheet                 [str]  CSV file containing the list of samples where the first column designates the file name of the
    	                                      normal sample, the second column for the file name of the matched tumor sample, example of the
    	                                      format for this file is in the testSamples directory
		-profile                       [str]  Configuration profile to use, each profile described in nextflow.config file
		                                      Available: preprocessing, germline, somatic

	Main Options:
		-bg                           [flag]  Runs the pipeline processes in the background, this option should be included if deploying
		                                      pipeline with real data set so processes will not be cut if user disconnects from deployment
		                                      environment
		-resume                       [flag]  Successfully completed tasks are cached so that if the pipeline stops prematurely the
		                                      previously completed tasks are skipped while maintaining their output
		--input_dir                    [str]  Directory that holds BAMs and associated index files
		                                      Default: ./input/preprocessedBams
		--output_dir                   [str]  Directory that will hold all output files from the somatic variant analysis
		                                      Default: ./output
		--email                        [str]  Email address to send workflow completion/stoppage notification
		--singularity_module           [str]  Indicates the name of the Singularity software module to be loaded for use in the pipeline,
		                                      this option is not needed if Singularity is natively installed on the deployment environment
		--mutect_ref_vcf_concatenated  [str]  Indicates whether or not the gnomAD allele frequency reference VCF used for MuTect2 processes has
		                                      been concatenated, this will be done in a process of the pipeline if it has not, this does not
		                                      need to be done for every separate run after the first
		                                      Available: yes, no
		                                      Default: yes
		--vep_ref_cached               [str]  Indicates whether or not the VEP reference files used for annotation have been downloaded/cached
		                                      locally, this will be done in a process of the pipeline if it has not, this does not need to be
		                                      done for every separate run after the first
		                                      Available: yes, no
		                                      Default: yes
		--cpus                         [int]  Globally set the number of cpus to be allocated for all processes
		                                      Available: 2, 4, 8, 16, etc.
		                                      Default: uniquly set for each process in ./nextflow.config to minimize resources needed
		--memory                       [str]  Globally set the amount of memory to be allocated for all processes, written as '##.GB' or '##.MB'
		                                      Available: 32.GB, 2400.MB, etc.
		                                      Default: uniquly set for each process in ./nextflow.config to minimize resources needed
		--help                        [flag]  Prints this message

	Toolbox Switches:
		--telomerehunter               [str]  Indicates whether or not to use this tool
		                                      Available: off, on
		                                      Default: on
		--conpair                      [str]  Indicates whether or not to use this tool
		                                      Available: off, on
		                                      Default: on
		--varscan                      [str]  Indicates whether or not to use this tool
		                                      Available: off, on
		                                      Default: on
		--mutect                       [str]  Indicates whether or not to use this tool
		                                      Available: off, on
		                                      Default: on
		--strelka                      [str]  Indicates whether or not to use this tool
		                                      Available: off, on
		                                      Default: on
		--caveman                      [str]  Indicates whether or not to use this tool
		                                      Available: off, on
		                                      Default: on
		--ascatngs                     [str]  Indicates whether or not to use this tool
		                                      Available: off, on
		                                      Default: on
		--controlfreec                 [str]  Indicates whether or not to use this tool
		                                      Available: off, on
		                                      Default: on
		--sclust                       [str]  Indicates whether or not to use this tool
		                                      Available: off, on
		                                      Default: on
		--manta                        [str]  Indicates whether or not to use this tool
		                                      Available: off, on
		                                      Default: on
		--svaba                        [str]  Indicates whether or not to use this tool
		                                      Available: off, on
		                                      Default: on
		--delly                        [str]  Indicates whether or not to use this tool
		                                      Available: off, on
		                                      Default: on
		--gridss DISABLED              [str]  Indicates whether or not to use this tool
		                                      Available: off, on
		                                      Default: off

	################################################

	""".stripIndent()
}

// #################################################### \\
// ~~~~~~~~~~~~~ PARAMETER CONFIGURATION ~~~~~~~~~~~~~~ \\

// Declare the defaults for all pipeline parameters
params.input_dir = "${workflow.projectDir}/input/preprocessedBams"
params.output_dir = "${workflow.projectDir}/output"
params.run_id = null
params.sample_sheet = null
params.mutect_ref_vcf_concatenated = "yes"
params.vep_ref_cached = "yes"
params.telomerehunter = "on"
params.conpair = "on"
params.varscan = "on"
params.mutect = "on"
params.strelka = "on"
params.caveman = "on"
params.ascatngs = "on"
params.controlfreec = "on"
params.sclust = "on"
params.manta = "on"
params.svaba = "on"
params.delly = "on"
params.gridss = "off"
params.ascatngs_ploidy = null
params.ascatngs_purity = null
params.cpus = null
params.memory = null
params.help = null

// Print help message if requested
if( params.help ) exit 0, helpMessage()

// Print erro message if user-defined input/output directories does not exist
if( !file(params.input_dir).exists() ) exit 1, "The user-specified input directory does not exist in filesystem."

// Print error messages if required parameters are not set
if( params.run_id == null ) exit 1, "The run command issued does not have the '--run_id' parameter set. Please set the '--run_id' parameter to a unique identifier for the run."

if( params.sample_sheet == null ) exit 1, "The run command issued does not have the '--sample_sheet' parameter set. Please set the '--sample_sheet' parameter to the path of the normal/tumor pair sample sheet CSV."

// Print preemptive error message if either ascatNGS ploidy or purity is set while the other is not
if( (params.ascatngs_ploidy && !params.ascatngs_purity) || (!params.ascatngs_ploidy && params.ascatngs_purity) ) exit 1, "User must define both ascatNGS ploidy and purity or leave both at default value"

// Set channels for reference files
Channel
	.fromPath( 'references/hg38/Homo_sapiens_assembly38.fasta' )
	.into{ reference_genome_fasta_forConpairPileup;
	       reference_genome_fasta_forVarscanSamtoolsMpileup;
	       reference_genome_fasta_forVarscanBamReadcount;
	       reference_genome_fasta_forVarscanBcftoolsNorm;
	       reference_genome_fasta_forMutectCalling;
	       reference_genome_fasta_forMutectFilter;
	       reference_genome_fasta_forMutectBcftools;
	       reference_genome_fasta_forCaveman;
	       reference_genome_fasta_forAscatNgs;
	       reference_genome_fasta_forControlFreecSamtoolsMpileup;
	       reference_genome_fasta_forControlFreecCalling;
	       reference_genome_fasta_forManta;
	       reference_genome_fasta_forStrelka;
	       reference_genome_fasta_forStrelkaBcftools;
	       reference_genome_fasta_forSvabaBcftools;
	       reference_genome_fasta_forDelly;
	       reference_genome_fasta_forGridssSetup;
	       reference_genome_fasta_forGridssPostprocessing;
	       reference_genome_fasta_forConsensusSnvMpileup;
	       reference_genome_fasta_forConsensusIndelMpileup;
	       reference_genome_fasta_forAnnotation }

Channel
	.fromPath( 'references/hg38/Homo_sapiens_assembly38.fasta.fai' )
	.into{ reference_genome_fasta_index_forAlleleCount;
		   reference_genome_fasta_index_forConpairPileup;
	       reference_genome_fasta_index_forVarscanSamtoolsMpileup;
	       reference_genome_fasta_index_forVarscanBamReadcount;
	       reference_genome_fasta_index_forVarscanBcftoolsNorm;
	       reference_genome_fasta_index_forMutectCalling;
	       reference_genome_fasta_index_forMutectFilter;
	       reference_genome_fasta_index_forMutectBcftools;
	       reference_genome_fasta_index_forCaveman;
	       reference_genome_fasta_index_forAscatNgs;
	       reference_genome_fasta_index_forControlFreecSamtoolsMpileup;
	       reference_genome_fasta_index_forControlFreecCalling;
	       reference_genome_fasta_index_forManta;
	       reference_genome_fasta_index_forStrelka;
	       reference_genome_fasta_index_forStrelkaBcftools;
	       reference_genome_fasta_index_forSvabaBcftools;
	       reference_genome_fasta_index_forDelly;
	       reference_genome_fasta_index_forGridssSetup;
	       reference_genome_fasta_index_forGridssPostprocessing;
	       reference_genome_fasta_index_forConsensusSnvMpileup;
	       reference_genome_fasta_index_forConsensusIndelMpileup;
	       reference_genome_fasta_index_forAnnotation }

Channel
	.fromPath( 'references/hg38/Homo_sapiens_assembly38.dict' )
	.into{ reference_genome_fasta_dict_forConpairPileup;
	       reference_genome_fasta_dict_forVarscanSamtoolsMpileup;
	       reference_genome_fasta_dict_forVarscanBamReadcount;
	       reference_genome_fasta_dict_forVarscanBcftoolsNorm;
	       reference_genome_fasta_dict_forMutectCalling;
	       reference_genome_fasta_dict_forMutectPileupGatherTumor;
	       reference_genome_fasta_dict_forMutectPileupGatherNormal;
	       reference_genome_fasta_dict_forMutectFilter;
	       reference_genome_fasta_dict_forMutectBcftools;
	       reference_genome_fasta_dict_forCaveman;
	       reference_genome_fasta_dict_forAscatNgs;
	       reference_genome_fasta_dict_forControlFreecSamtoolsMpileup;
	       reference_genome_fasta_dict_forControlFreecCalling;
	       reference_genome_fasta_dict_forManta;
	       reference_genome_fasta_dict_forStrelka;
	       reference_genome_fasta_dict_forStrelkaBcftools;
	       reference_genome_fasta_dict_forSvaba;
	       reference_genome_fasta_dict_forSvabaBcftools;
	       reference_genome_fasta_dict_forDelly;
	       reference_genome_fasta_dict_forGridssPostprocessing;
	       reference_genome_fasta_dict_forConsensusSnvMpileup;
	       reference_genome_fasta_dict_forConsensusIndelMpileup;
	       reference_genome_fasta_dict_forAnnotation }

Channel
	.fromPath( 'references/hg38/wgs_calling_regions.hg38.bed' )
	.into{ gatk_bundle_wgs_bed_forVarscanSamtoolsMpileup;
	       gatk_bundle_wgs_bed_forMutectCalling;
	       gatk_bundle_wgs_bed_forMutectPileup;
	       gatk_bundle_wgs_bed_forControlFreecSamtoolsMpileup;
	       gatk_bundle_wgs_bed_forManta;
	       gatk_bundle_wgs_bed_forStrelka;
	       gatk_bundle_wgs_bed_forSvaba }

Channel
	.fromPath( 'references/hg38/wgs_calling_regions_blacklist.0based.hg38.bed' )
	.set{ gatk_bundle_wgs_bed_blacklist_0based_forDelly }

Channel
	.fromPath( 'references/hg38/wgs_calling_regions_blacklist.1based.hg38.bed' )
	.into{ gatk_bundle_wgs_bed_blacklist_1based_forCaveman;
	       gatk_bundle_wgs_bed_blacklist_1based_forGridssSetup }

Channel
	.fromList( ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6',
	            'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
	            'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18',
	            'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY',] )
	.into{ chromosome_list_forVarscanSamtoolsMpileup;
	       chromosome_list_forMutectCalling;
	       chromosome_list_forMutectPileup;
	       chromosome_list_forCavemanSplit;
	       chromosome_list_forCavemanMstep;
	       chromosome_list_forCavemanEstep;
	       chromosome_list_forControlFreecSamtoolsMpileup;
	       chromosome_list_forControlFreecMerge;
	       chromosome_list_forSclustBamprocess }

Channel
	.fromPath( 'references/hg38/sex_identification_loci.chrY.hg38.txt' )
	.set{ sex_identification_loci }

Channel
	.fromPath( 'references/hg38/cytoband_autosome_sex_chroms.hg38.bed' )
	.set{ cytoband_bed }

Channel
	.fromPath( 'references/hg38/1000g_pon.hg38.vcf.gz' )
	.set{ panel_of_normals_1000G }

Channel
	.fromPath( 'references/hg38/1000g_pon.hg38.vcf.gz.tbi' )
	.set{ panel_of_normals_1000G_index }	

Channel
	.fromPath( 'references/hg38/af-only-gnomad.chr1-9.hg38.vcf.gz' )
	.set{ gnomad_ref_vcf_chromosomes1_9 }

Channel
	.fromPath( 'references/hg38/af-only-gnomad.chr1-9.hg38.vcf.gz.tbi' )
	.set{ gnomad_ref_vcf_chromosomes1_9_index }

Channel
	.fromPath( 'references/hg38/af-only-gnomad.chr10-22.hg38.vcf.gz' )
	.set{ gnomad_ref_vcf_chromosomes10_22 }

Channel
	.fromPath( 'references/hg38/af-only-gnomad.chr10-22.hg38.vcf.gz.tbi' )
	.set{ gnomad_ref_vcf_chromosomes10_22_index }

Channel
	.fromPath( 'references/hg38/af-only-gnomad.chrXYM-alts.hg38.vcf.gz' )
	.set{ gnomad_ref_vcf_chromosomesXYM_alts }

Channel
	.fromPath( 'references/hg38/af-only-gnomad.chrXYM-alts.hg38.vcf.gz.tbi' )
	.set{ gnomad_ref_vcf_chromosomesXYM_alts_index }

if( params.mutect == "on" && params.mutect_ref_vcf_concatenated == "yes" ) {
	Channel
		.fromPath( 'references/hg38/af-only-gnomad.hg38.vcf.gz', checkIfExists: true )
		.ifEmpty{ error "The run command issued has the '--mutect_ref_vcf_concatenated' parameter set to 'yes', however the file does not exist. Please set the '--mutect_ref_vcf_concatenated' parameter to 'no' and resubmit the run command. For more information, check the README or issue the command 'nextflow run somatic.nf --help'"}
		.set{ mutect_gnomad_ref_vcf_preBuilt }

	Channel
		.fromPath( 'references/hg38/af-only-gnomad.hg38.vcf.gz.tbi', checkIfExists: true )
		.ifEmpty{ error "The '--mutect_ref_vcf_concatenated' parameter set to 'yes', however the index file does not exist for the reference VCF. Please set the '--mutect_ref_vcf_concatenated' parameter to 'no' and resubmit the run command. Alternatively, use Tabix to index the reference VCF."}
		.set{ mutect_gnomad_ref_vcf_index_preBuilt }
}

Channel
	.fromPath( 'references/hg38/small_exac_common_3.hg38.vcf.gz' )
	.set{ exac_common_sites_ref_vcf }

Channel
	.fromPath( 'references/hg38/small_exac_common_3.hg38.vcf.gz.tbi' )
	.set{ exac_common_sites_ref_vcf_index }

Channel
	.fromPath( 'references/hg38/SnpGcCorrections.hg38.tsv' )
	.set{ snp_gc_corrections }

Channel
	.fromPath( 'references/hg38/Homo_sapiens_assembly38_autosome_sex_chroms', type: 'dir' )
	.set{ autosome_sex_chromosome_fasta_dir }

Channel
	.fromPath( 'references/hg38/Homo_sapiens_assembly38_autosome_sex_chrom_sizes.txt' )
	.set{ autosome_sex_chromosome_sizes }

Channel
	.fromPath( 'references/hg38/common_all_20180418.vcf.gz' )
	.set{ common_dbsnp_ref_vcf }

Channel
	.fromPath( 'references/hg38/common_all_20180418.vcf.gz.tbi' )
	.set{ common_dbsnp_ref_vcf_index }

Channel
	.fromPath( 'references/hg38/mappability_track_100m2.hg38.zip' )
	.set{ mappability_track_zip }

Channel
	.fromPath( ['references/hg38/Homo_sapiens_assembly38.fasta', 'references/hg38/Homo_sapiens_assembly38.fasta.fai',
	            'references/hg38/Homo_sapiens_assembly38.fasta.64.alt', 'references/hg38/Homo_sapiens_assembly38.fasta.64.amb',
	            'references/hg38/Homo_sapiens_assembly38.fasta.64.ann', 'references/hg38/Homo_sapiens_assembly38.fasta.64.bwt',
	            'references/hg38/Homo_sapiens_assembly38.fasta.64.pac', 'references/hg38/Homo_sapiens_assembly38.fasta.64.sa'] )
	.set{ bwa_ref_genome_files }

Channel
	.fromPath( 'references/hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz' )
	.set{ dbsnp_known_indel_ref_vcf }

Channel
	.fromPath( 'references/hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi' )
	.set{ dbsnp_known_indel_ref_vcf_index }

Channel
	.fromPath( 'references/hg38/pon_single_breakend.hg38.bed' )
	.set{ pon_single_breakend_bed }

Channel
	.fromPath( 'references/hg38/pon_breakpoint.hg38.bedpe' )
	.set{ pon_breakpoint_bedpe }

Channel
	.fromPath( 'references/hg38/known_fusion_pairs_v3.hg38.bedpe' )
	.set{ known_fusion_pairs_bedpe }

Channel
	.fromPath( 'references/hg38/centromeric_repeats.hg38.bed.gz' )
	.set{ centromeric_repeats_bed }

Channel
	.fromPath( 'references/hg38/centromeric_repeats.hg38.bed.gz.tbi' )
	.set{ centromeric_repeats_bed_index }

Channel
	.fromPath( 'references/hg38/simple_repeats.hg38.bed.gz' )
	.set{ simple_repeats_bed }

Channel
	.fromPath( 'references/hg38/simple_repeats.hg38.bed.gz.tbi' )
	.set{ simple_repeats_bed_index }

Channel
	.fromPath( 'references/hg38/dbsnp138.hg38.bed.gz' )
	.set{ dbsnp_bed }

Channel
	.fromPath( 'references/hg38/dbsnp138.hg38.bed.gz.tbi' )
	.set{ dbsnp_bed_index }

Channel
	.fromPath( 'references/hg38/unmatched_normal.cpbi.hg38.bed.gz' )
	.set{ unmatched_normal_bed }

Channel
	.fromPath( 'references/hg38/unmatched_normal.cpbi.hg38.bed.gz.tbi' )
	.set{ unmatched_normal_bed_index }

if( params.vep_ref_cached == "yes" ) {
	Channel
		.fromPath( 'references/hg38/homo_sapiens_vep_101_GRCh38/', type: 'dir', checkIfExists: true )
		.ifEmpty{ error "The run command issued has the '--vep_ref_cached' parameter set to 'yes', however the directory does not exist. Please set the '--vep_ref_cached' parameter to 'no' and resubmit the run command. For more information, check the README or issue the command 'nextflow run somatic.nf --help'"}
		.set{ vep_ref_dir_preDownloaded }
}

// #################################################### \\
// ~~~~~~~~~~~~~~~~ PIPELINE PROCESSES ~~~~~~~~~~~~~~~~ \\

log.info ''
log.info '##### Myeloma Genome Project 1000 Pipeline #####'
log.info '################################################'
log.info '~~~~~~~~~~~ SOMATIC VARIANT ANALYSIS ~~~~~~~~~~~'
log.info '################################################'
log.info ''
log.info "~~~ Launch Time ~~~		${workflowTimestamp}"
log.info ''
log.info "~~~ Input Directory ~~~ 	${params.input_dir}"
log.info ''
log.info "~~~ Output Directory ~~~ 	${params.output_dir}"
log.info ''
log.info "~~~ Run Report File ~~~ 	nextflow_report.${params.run_id}.html"
log.info ''
log.info '################################################'
log.info ''

// Read user provided sample sheet to set the Tumor/Normal sample pairs
Channel
	.fromPath( params.sample_sheet )
	.splitCsv( header:true )
	.map{ row -> tumor_bam = "${row.tumor}"
				 tumor_bam_index = "${row.tumor}".replaceFirst(/\.bam$/, "")
	             normal_bam = "${row.normal}"
	             normal_bam_index = "${row.normal}".replaceFirst(/\.bam$/, "")
	             return[ file("${params.input_dir}/${tumor_bam}"),
	             		 file("${params.input_dir}/${tumor_bam_index}*.bai"),
	             		 file("${params.input_dir}/${normal_bam}"), 
	             		 file("${params.input_dir}/${normal_bam_index}*.bai") ] }
	.into{ tumor_normal_pair_forAlleleCount;
		   tumor_normal_pair_forTelomereHunter;
		   tumor_normal_pair_forConpairPileup;
	       tumor_normal_pair_forVarscanSamtoolsMpileup;
	       tumor_normal_pair_forMutectCalling;
	       tumor_normal_pair_forMutectPileup;
	       tumor_normal_pair_forControlFreecSamtoolsMpileup;
	       tumor_normal_pair_forSclustBamprocess;
	       tumor_normal_pair_forManta;
	       tumor_normal_pair_forSvaba;
	       tumor_normal_pair_forDelly;
	       tumor_normal_pair_forGridssPreprocess }

// Combine reference FASTA index and sex identification loci files into one channel for use in alleleCount process
reference_genome_fasta_index_forAlleleCount.combine( sex_identification_loci )
	.set{ ref_index_and_sex_ident_loci }

// alleleCount ~ determine the sex of each sample to use in downstream analyses
process identifySampleSex_allelecount {
	publishDir "${params.output_dir}/somatic/sexOfSamples", mode: 'copy', pattern: '*.{txt}'
	tag "T=${tumor_id} N=${normal_id}"

	input:
	tuple path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(reference_genome_fasta_index_forAlleleCount), path(sex_identification_loci) from tumor_normal_pair_forAlleleCount.combine(ref_index_and_sex_ident_loci)

	output:
	tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index) into bams_forVarscanBamReadcount
	tuple val(tumor_normal_sample_id), path(sample_sex) into sex_of_sample_forControlFreecCalling
	tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(sample_sex) into bams_and_sex_of_sample_forAscatNgs
	tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index) into bams_forConsensusSnvMpileup
	tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index) into bams_forConsensusIndelMpileup

	script:
	tumor_id = "${tumor_bam.baseName}".replaceFirst(/\..*$/, "")
	normal_id = "${normal_bam.baseName}".replaceFirst(/\..*$/, "")
	tumor_normal_sample_id = "${tumor_id}_vs_${normal_id}"
	sex_loci_allele_counts = "${tumor_normal_sample_id}.sexloci.txt"
	sample_sex = "${tumor_normal_sample_id}.sexident.txt"
	"""
	alleleCounter \
	--loci-file "${sex_identification_loci}" \
	--hts-file "${normal_bam}" \
	--ref-file "${reference_genome_fasta_index_forAlleleCount}" \
	--output-file "${sex_loci_allele_counts}"

	sample_sex_determinator.sh "${sex_loci_allele_counts}" > "${sample_sex}"
	"""
}

// TelomereHunter ~ estimate telomere content and composition
process telomereEstimation_telomerehunter {
	publishDir "${params.output_dir}/somatic/telomereHunter", mode: 'copy'
	tag "${tumor_normal_sample_id}"

	input:
	tuple path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(cytoband_bed) from tumor_normal_pair_forTelomereHunter.combine(cytoband_bed)

	output:
	path "${tumor_normal_sample_id}/*.tsv"
	path "${tumor_normal_sample_id}/*.png"
	path "${tumor_normal_sample_id}/control_TelomerCnt_${tumor_normal_sample_id}/*.tsv"
	path "${tumor_normal_sample_id}/control_TelomerCnt_${tumor_normal_sample_id}/TVRs"
	path "${tumor_normal_sample_id}/tumor_TelomerCnt_${tumor_normal_sample_id}/*.tsv"
	path "${tumor_normal_sample_id}/tumor_TelomerCnt_${tumor_normal_sample_id}/TVRs"
	path "${tumor_normal_sample_id}/plots"

	when:
	params.telomerehunter == "on"

	script:
	tumor_id = "${tumor_bam.baseName}".replaceFirst(/\..*$/, "")
	normal_id = "${normal_bam.baseName}".replaceFirst(/\..*$/, "")
	tumor_normal_sample_id = "${tumor_id}_vs_${normal_id}"
	"""
	telomerehunter \
	--inputBamTumor "${tumor_bam}" \
	--inputBamControl "${normal_bam}" \
	--outPath . \
	--pid "${tumor_normal_sample_id}" \
	--bandingFile "${cytoband_bed}" \
	--parallel \
	--plotFileFormat png
	"""
}

// ~~~~~~~~~~~~~~~~ Conpair ~~~~~~~~~~~~~~~~ \\
// START

// Combine all reference FASTA files into one channel for use in Conpair Pileup process
reference_genome_fasta_forConpairPileup.combine( reference_genome_fasta_index_forConpairPileup )
	.combine( reference_genome_fasta_dict_forConpairPileup )
	.set{ reference_genome_bundle_forConpairPileup }

// Conpair run_gatk_pileup_for_sample ~ generate GATK pileups the tumor and normal BAMs separately
process bamPileupForConpair_conpair {
	publishDir "${params.output_dir}/somatic/conpair/pileups", mode: 'symlink'
	tag "T=${tumor_id} N=${normal_id}"

	input:
	tuple path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(reference_genome_fasta_forConpairPileup), path(reference_genome_fasta_index_forConpairPileup), path(reference_genome_fasta_dict_forConpairPileup) from tumor_normal_pair_forConpairPileup.combine(reference_genome_bundle_forConpairPileup)

	output:
	tuple val(tumor_id), val(normal_id) into sample_ids_forConpair
	tuple path(tumor_pileup), path(normal_pileup) into bam_pileups_forConpair

	when:
	params.conpair == "on"

	script:
	tumor_id = "${tumor_bam.baseName}".replaceFirst(/\..*$/, "")
	normal_id = "${normal_bam.baseName}".replaceFirst(/\..*$/, "")
	tumor_normal_sample_id = "${tumor_id}_vs_${normal_id}"
	tumor_pileup = "${tumor_id}.pileup"
	normal_pileup = "${normal_id}.pileup"
	hg38_ref_genome_markers = "/data/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover"
	"""
	\${CONPAIR_DIR}/scripts/run_gatk_pileup_for_sample.py \
	--xmx_jav "${task.memory.toGiga()}g" \
	--bam "${tumor_bam}" \
	--outfile "${tumor_pileup}" \
	--reference "${reference_genome_fasta_forConpairPileup}" \
	--markers \${CONPAIR_DIR}"${hg38_ref_genome_markers}.bed"

	\${CONPAIR_DIR}/scripts/run_gatk_pileup_for_sample.py \
	--xmx_jav "${task.memory.toGiga()}g" \
	--bam "${normal_bam}" \
	--outfile "${normal_pileup}" \
	--reference "${reference_genome_fasta_forConpairPileup}" \
	--markers \${CONPAIR_DIR}"${hg38_ref_genome_markers}.bed"
	"""
}

// Conpair verify_concordance / estimate_tumor_normal_contamination ~ concordance and contamination estimator for tumor–normal pileups
process concordanceAndContaminationEstimation_conpair {
	publishDir "${params.output_dir}/somatic/conpair", mode: 'copy'
	tag "T=${tumor_id} N=${normal_id}"

	input:
	tuple val(tumor_id), val(normal_id) from sample_ids_forConpair
	tuple path(tumor_pileup), path(normal_pileup) from bam_pileups_forConpair

	output:
	path concordance_file
	path contamination_file

	when:
	params.conpair == "on"
	
	script:
	concordance_file = "${tumor_id}_vs_${normal_id}.concordance.txt"
	contamination_file = "${tumor_id}_vs_${normal_id}.contamination.txt"
	hg38_ref_genome_markers = "/data/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover"
	"""
	\${CONPAIR_DIR}/scripts/verify_concordance.py \
	--min_cov 10 \
	--min_mapping_quality 10 \
	--min_base_quality 20 \
	--tumor_pileup "${tumor_pileup}" \
	--normal_pileup "${normal_pileup}" \
	--outfile "${concordance_file}" \
	--markers \${CONPAIR_DIR}"${hg38_ref_genome_markers}.txt"

	\${CONPAIR_DIR}/scripts/estimate_tumor_normal_contamination.py \
	--grid 0.01 \
	--min_mapping_quality 10 \
	--tumor_pileup "${tumor_pileup}" \
	--normal_pileup "${normal_pileup}" \
	--outfile "${contamination_file}" \
	--markers \${CONPAIR_DIR}"${hg38_ref_genome_markers}.txt"
	"""
}

// END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \\


// ~~~~~~~~~~~~~~~~ VarScan2 ~~~~~~~~~~~~~~~~ \\
// START

// Combine all reference FASTA files and WGS BED file into one channel for use in VarScan / SAMtools mpileup
reference_genome_fasta_forVarscanSamtoolsMpileup.combine( reference_genome_fasta_index_forVarscanSamtoolsMpileup )
	.combine( reference_genome_fasta_dict_forVarscanSamtoolsMpileup )
	.set{ reference_genome_bundle_forVarscanSamtoolsMpileup }

reference_genome_bundle_forVarscanSamtoolsMpileup.combine( gatk_bundle_wgs_bed_forVarscanSamtoolsMpileup )
	.set{ reference_genome_bundle_and_bed_forVarscanSamtoolsMpileup }

// VarScan somatic / SAMtools mpileup ~ heuristic/statistic approach to call SNV and indel variants
process snvAndIndelCalling_varscan {
	publishDir "${params.output_dir}/somatic/varscan/rawVcfs", mode: 'symlink'
	tag "C=${chromosome} T=${tumor_id} N=${normal_id}"

	input:
	tuple path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(reference_genome_fasta_forVarscanSamtoolsMpileup), path(reference_genome_fasta_index_forVarscanSamtoolsMpileup), path(reference_genome_fasta_dict_forVarscanSamtoolsMpileup), path(gatk_bundle_wgs_bed_forVarscanSamtoolsMpileup) from tumor_normal_pair_forVarscanSamtoolsMpileup.combine(reference_genome_bundle_and_bed_forVarscanSamtoolsMpileup)
	each chromosome from chromosome_list_forVarscanSamtoolsMpileup

	output:
	tuple val(tumor_normal_sample_id), path(raw_per_chromosome_snv_vcf), path(raw_per_chromosome_snv_vcf_index), path(raw_per_chromosome_indel_vcf), path(raw_per_chromosome_indel_vcf_index) into raw_per_chromosome_vcfs_forVarscanBcftools

	when:
	params.varscan == "on"

	script:
	tumor_id = "${tumor_bam.baseName}".replaceFirst(/\..*$/, "")
	normal_id = "${normal_bam.baseName}".replaceFirst(/\..*$/, "")
	tumor_normal_sample_id = "${tumor_id}_vs_${normal_id}"
	raw_per_chromosome_snv_vcf = "${tumor_normal_sample_id}.${chromosome}.snv.vcf.gz"
	raw_per_chromosome_snv_vcf_index = "${raw_per_chromosome_snv_vcf}.tbi"
	raw_per_chromosome_indel_vcf = "${tumor_normal_sample_id}.${chromosome}.indel.vcf.gz"
	raw_per_chromosome_indel_vcf_index = "${raw_per_chromosome_indel_vcf}.tbi"
	"""
	samtools mpileup \
	--no-BAQ \
	--min-MQ 1 \
	--positions "${gatk_bundle_wgs_bed_forVarscanSamtoolsMpileup}" \
	--region "${chromosome}" \
	--fasta-ref "${reference_genome_fasta_forVarscanSamtoolsMpileup}" \
	"${normal_bam}" "${tumor_bam}" \
	| \
	java -jar \${VARSCAN} somatic \
	--mpileup 1 \
	--min-coverage-normal 8 \
	--min-coverage-tumor 6 \
	--min-var-freq 0.10 \
	--min-freq-for-hom 0.75 \
	--normal-purity 1.00 \
	--tumor-purity 1.00 \
	--p-value 0.99 \
	--somatic-p-value 0.05 \
	--strand-filter 0 \
	--output-vcf \
	--output-snp "${tumor_normal_sample_id}.${chromosome}.snv" \
	--output-indel "${tumor_normal_sample_id}.${chromosome}.indel"

	bgzip < "${tumor_normal_sample_id}.${chromosome}.snv.vcf" > "${raw_per_chromosome_snv_vcf}"
	tabix "${raw_per_chromosome_snv_vcf}"

	bgzip < "${tumor_normal_sample_id}.${chromosome}.indel.vcf" > "${raw_per_chromosome_indel_vcf}"
	tabix "${raw_per_chromosome_indel_vcf}"
	"""
}

// BCFtools concat ~ concatenate all VarScan SNV/indel per chromosome VCFs
process concatenateVarscanPerChromosomeVcfs_bcftools {
	publishDir "${params.output_dir}/somatic/varscan/rawVcfs", mode: 'symlink'
	tag "${tumor_normal_sample_id}"
	
	input:
	tuple val(tumor_normal_sample_id), path(raw_per_chromosome_snv_vcf), path(raw_per_chromosome_snv_vcf_index), path(raw_per_chromosome_indel_vcf), path(raw_per_chromosome_indel_vcf_index) from raw_per_chromosome_vcfs_forVarscanBcftools.groupTuple()

	output:
	tuple val(tumor_normal_sample_id), path(raw_snv_vcf), path(raw_indel_vcf) into raw_vcfs_forVarscanHcFilter

	when:
	params.varscan == "on"

	script:
	raw_snv_vcf = "${tumor_normal_sample_id}.snv.vcf.gz"
	raw_indel_vcf = "${tumor_normal_sample_id}.indel.vcf.gz"
	"""
	bcftools concat \
	--threads ${task.cpus} \
	--output-type z \
	--output "${raw_snv_vcf}" \
	"${tumor_normal_sample_id}.chr1.snv.vcf.gz" \
	"${tumor_normal_sample_id}.chr2.snv.vcf.gz" \
	"${tumor_normal_sample_id}.chr3.snv.vcf.gz" \
	"${tumor_normal_sample_id}.chr4.snv.vcf.gz" \
	"${tumor_normal_sample_id}.chr5.snv.vcf.gz" \
	"${tumor_normal_sample_id}.chr6.snv.vcf.gz" \
	"${tumor_normal_sample_id}.chr7.snv.vcf.gz" \
	"${tumor_normal_sample_id}.chr8.snv.vcf.gz" \
	"${tumor_normal_sample_id}.chr9.snv.vcf.gz" \
	"${tumor_normal_sample_id}.chr10.snv.vcf.gz" \
	"${tumor_normal_sample_id}.chr11.snv.vcf.gz" \
	"${tumor_normal_sample_id}.chr12.snv.vcf.gz" \
	"${tumor_normal_sample_id}.chr13.snv.vcf.gz" \
	"${tumor_normal_sample_id}.chr14.snv.vcf.gz" \
	"${tumor_normal_sample_id}.chr15.snv.vcf.gz" \
	"${tumor_normal_sample_id}.chr16.snv.vcf.gz" \
	"${tumor_normal_sample_id}.chr17.snv.vcf.gz" \
	"${tumor_normal_sample_id}.chr18.snv.vcf.gz" \
	"${tumor_normal_sample_id}.chr19.snv.vcf.gz" \
	"${tumor_normal_sample_id}.chr20.snv.vcf.gz" \
	"${tumor_normal_sample_id}.chr21.snv.vcf.gz" \
	"${tumor_normal_sample_id}.chr22.snv.vcf.gz" \
	"${tumor_normal_sample_id}.chrX.snv.vcf.gz" \
	"${tumor_normal_sample_id}.chrY.snv.vcf.gz"

	bcftools concat \
	--threads ${task.cpus} \
	--output-type z \
	--output "${raw_indel_vcf}" \
	"${tumor_normal_sample_id}.chr1.indel.vcf.gz" \
	"${tumor_normal_sample_id}.chr2.indel.vcf.gz" \
	"${tumor_normal_sample_id}.chr3.indel.vcf.gz" \
	"${tumor_normal_sample_id}.chr4.indel.vcf.gz" \
	"${tumor_normal_sample_id}.chr5.indel.vcf.gz" \
	"${tumor_normal_sample_id}.chr6.indel.vcf.gz" \
	"${tumor_normal_sample_id}.chr7.indel.vcf.gz" \
	"${tumor_normal_sample_id}.chr8.indel.vcf.gz" \
	"${tumor_normal_sample_id}.chr9.indel.vcf.gz" \
	"${tumor_normal_sample_id}.chr10.indel.vcf.gz" \
	"${tumor_normal_sample_id}.chr11.indel.vcf.gz" \
	"${tumor_normal_sample_id}.chr12.indel.vcf.gz" \
	"${tumor_normal_sample_id}.chr13.indel.vcf.gz" \
	"${tumor_normal_sample_id}.chr14.indel.vcf.gz" \
	"${tumor_normal_sample_id}.chr15.indel.vcf.gz" \
	"${tumor_normal_sample_id}.chr16.indel.vcf.gz" \
	"${tumor_normal_sample_id}.chr17.indel.vcf.gz" \
	"${tumor_normal_sample_id}.chr18.indel.vcf.gz" \
	"${tumor_normal_sample_id}.chr19.indel.vcf.gz" \
	"${tumor_normal_sample_id}.chr20.indel.vcf.gz" \
	"${tumor_normal_sample_id}.chr21.indel.vcf.gz" \
	"${tumor_normal_sample_id}.chr22.indel.vcf.gz" \
	"${tumor_normal_sample_id}.chrX.indel.vcf.gz" \
	"${tumor_normal_sample_id}.chrY.indel.vcf.gz"
	"""
}

// VarScan processSomatic ~ filter the called SNVs and indels for confidence and somatic status assignment
process filterRawSnvAndIndels_varscan {
	publishDir "${params.output_dir}/somatic/varscan/filteredVcfs", mode: 'symlink'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(raw_snv_vcf), path(raw_indel_vcf) from raw_vcfs_forVarscanHcFilter

	output:
	tuple val(tumor_normal_sample_id), path(high_confidence_snv_vcf), path(high_confidence_snv_vcf_index), path(high_confidence_indel_vcf), path(high_confidence_indel_vcf_index) into high_confidence_vcfs_forVarscanBamReadcount, high_confidence_vcfs_forVarscanFpFilter

	when:
	params.varscan == "on"

	script:
	high_confidence_snv_vcf = "${raw_snv_vcf}".replaceFirst(/\.vcf\.gz/, ".somatic.hc.vcf.gz")
	high_confidence_snv_vcf_index = "${high_confidence_snv_vcf}.tbi"
	high_confidence_indel_vcf = "${raw_indel_vcf}".replaceFirst(/\.vcf\.gz/, ".somatic.hc.vcf.gz")
	high_confidence_indel_vcf_index = "${high_confidence_indel_vcf}.tbi"
	"""
	zcat "${raw_snv_vcf}" \
	| \
	java -jar \${VARSCAN} processSomatic \
	"${tumor_normal_sample_id}.snv" \
	--min-tumor-freq 0.10 \
	--max-normal-freq 0.05 \
	--p-value 0.07

	bgzip < "${tumor_normal_sample_id}.snv.Somatic.hc" > "${high_confidence_snv_vcf}"
	tabix "${high_confidence_snv_vcf}"

	zcat "${raw_indel_vcf}" \
	| \
	java -jar \${VARSCAN} processSomatic \
	"${tumor_normal_sample_id}.indel" \
	--min-tumor-freq 0.10 \
	--max-normal-freq 0.05 \
	--p-value 0.07

	bgzip < "${tumor_normal_sample_id}.indel.Somatic.hc" > "${high_confidence_indel_vcf}"
	tabix "${high_confidence_indel_vcf}"
	"""
}

// Combine all needed reference FASTA files into one channel for use in bam-readcount process
reference_genome_fasta_forVarscanBamReadcount.combine( reference_genome_fasta_index_forVarscanBamReadcount )
	.combine( reference_genome_fasta_dict_forVarscanBamReadcount )
	.set{ reference_genome_bundle_forVarscanBamReadcount }

// bam-readcount / BCFtools concat ~ generate metrics at single nucleotide positions for filtering out false positive calls
process bamReadcountForVarscanFpFilter_bamreadcount {
	publishDir "${params.output_dir}/somatic/varscan/filteredVcfs", mode: 'symlink'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(high_confidence_snv_vcf), path(high_confidence_snv_vcf_index), path(high_confidence_indel_vcf), path(high_confidence_indel_vcf_index), path(reference_genome_fasta_forVarscanBamReadcount), path(reference_genome_fasta_index_forVarscanBamReadcount), path(reference_genome_fasta_dict_forVarscanBamReadcount) from bams_forVarscanBamReadcount.join(high_confidence_vcfs_forVarscanBamReadcount).combine(reference_genome_bundle_forVarscanBamReadcount)

	output:
	tuple val(tumor_normal_sample_id), path(snv_readcount_file), path(indel_readcount_file) into readcount_forVarscanFpFilter

	when:
	params.varscan == "on"

	script:
	snv_readcount_file = "${tumor_normal_sample_id}_bam_readcount_snv.tsv"
	indel_readcount_file = "${tumor_normal_sample_id}_bam_readcount_indel.tsv"
	"""
	bcftools concat \
	--threads ${task.cpus} \
	--allow-overlaps \
	--output-type z \
	--output "${tumor_normal_sample_id}.somatic.hc.vcf.gz" \
	"${high_confidence_snv_vcf}" "${high_confidence_indel_vcf}"

	tabix "${tumor_normal_sample_id}.somatic.hc.vcf.gz"

	bam_readcount_helper.py \
	"${tumor_normal_sample_id}.somatic.hc.vcf.gz" \
	TUMOR \
	"${reference_genome_fasta_forVarscanBamReadcount}" \
	"${tumor_bam}" \
	.

	mv TUMOR_bam_readcount_snv.tsv "${snv_readcount_file}"
	mv TUMOR_bam_readcount_indel.tsv "${indel_readcount_file}"
	"""
}

// VarScan fpfilter ~ filter out additional false positive variants
process falsePositivefilterSnvAndIndels_varscan {
	publishDir "${params.output_dir}/somatic/varscan/filteredVcfs", mode: 'symlink'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(high_confidence_snv_vcf), path(high_confidence_snv_vcf_index), path(high_confidence_indel_vcf), path(high_confidence_indel_vcf_index), path(snv_readcount_file), path(indel_readcount_file) from high_confidence_vcfs_forVarscanFpFilter.join(readcount_forVarscanFpFilter)
	
	output:
	tuple val(tumor_normal_sample_id), path(fp_filtered_snv_vcf), path(fp_filtered_indel_vcf) into filtered_vcfs_forVarscanBcftools

	when:
	params.varscan == "on"

	script:
	unzipped_hc_snv_vcf = "${high_confidence_snv_vcf}".replaceFirst(/\.gz/, "")
	unzipped_hc_indel_vcf = "${high_confidence_indel_vcf}".replaceFirst(/\.gz/, "")
	fp_filtered_snv_vcf = "${unzipped_hc_snv_vcf}".replaceFirst(/\.hc\.vcf/, ".filtered.vcf")
	fp_filtered_indel_vcf = "${unzipped_hc_indel_vcf}".replaceFirst(/\.hc\.vcf/, ".filtered.vcf")
	"""
	gunzip -f "${high_confidence_snv_vcf}"

	java -jar \$VARSCAN fpfilter \
	"${unzipped_hc_snv_vcf}" \
	"${snv_readcount_file}" \
	--filtered-file "${tumor_normal_sample_id}.snv.failed.vcf" \
	--min-var-count 2 \
	--min-var-freq 0.01 \
	--min-ref-basequal 25 \
	--min-var-basequal 25 \
	--output-file "${fp_filtered_snv_vcf}"

	gunzip -f "${high_confidence_indel_vcf}"

	java -jar \$VARSCAN fpfilter \
	"${unzipped_hc_indel_vcf}" \
	"${indel_readcount_file}" \
	--filtered-file "${tumor_normal_sample_id}.indel.failed.vcf" \
	--min-var-count 2 \
	--min-var-freq 0.01 \
	--min-ref-basequal 25 \
	--min-var-basequal 25 \
	--output-file "${fp_filtered_indel_vcf}"
	"""
}

// Combine all needed reference FASTA files into one channel for use in VarScan / BCFtools Norm process
reference_genome_fasta_forVarscanBcftoolsNorm.combine( reference_genome_fasta_index_forVarscanBcftoolsNorm )
	.combine( reference_genome_fasta_dict_forVarscanBcftoolsNorm )
	.set{ reference_genome_bundle_forVarscanBcftoolsNorm }

// BCFtools norm ~ split multiallelic sites into multiple rows then left-align and normalize indels
process splitMultiallelicAndLeftNormalizeVarscanVcf_bcftools {
	publishDir "${params.output_dir}/somatic/varscan", mode: 'copy', pattern: '*.{vcf.gz,tbi,txt}'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(fp_filtered_snv_vcf), path(fp_filtered_indel_vcf), path(reference_genome_fasta_forVarscanBcftoolsNorm), path(reference_genome_fasta_index_forVarscanBcftoolsNorm), path(reference_genome_fasta_dict_forVarscanBcftoolsNorm) from filtered_vcfs_forVarscanBcftools.combine(reference_genome_bundle_forVarscanBcftoolsNorm)

	output:
	tuple val(tumor_normal_sample_id), path(final_varscan_snv_vcf), path(final_varscan_snv_vcf_index) into final_varscan_snv_vcf_forConsensus
	tuple val(tumor_normal_sample_id), path(final_varscan_indel_vcf), path(final_varscan_indel_vcf_index) into final_varscan_indel_vcf_forConsensus
	path varscan_snv_multiallelics_stats
	path varscan_indel_multiallelics_stats
	path varscan_indel_realign_normalize_stats

	when:
	params.varscan == "on"

	script:
	final_varscan_snv_vcf = "${tumor_normal_sample_id}.varscan.somatic.snv.vcf.gz"
	final_varscan_snv_vcf_index ="${final_varscan_snv_vcf}.tbi"
	final_varscan_indel_vcf = "${tumor_normal_sample_id}.varscan.somatic.indel.vcf.gz"
	final_varscan_indel_vcf_index = "${final_varscan_indel_vcf}.tbi"
	varscan_snv_multiallelics_stats = "${tumor_normal_sample_id}.varscan.snv.multiallelicsstats.txt"
	varscan_indel_multiallelics_stats = "${tumor_normal_sample_id}.varscan.indel.multiallelicsstats.txt"
	varscan_indel_realign_normalize_stats = "${tumor_normal_sample_id}.varscan.indel.realignnormalizestats.txt"
	"""
	bgzip --stdout < "${fp_filtered_snv_vcf}" \
	| \
	bcftools norm \
	--threads ${task.cpus} \
	--multiallelics -snps \
	--output-type z \
	--output "${final_varscan_snv_vcf}" \
	- 2>"${varscan_snv_multiallelics_stats}"

	tabix "${final_varscan_snv_vcf}"
	
	bgzip --stdout < "${fp_filtered_indel_vcf}" \
	| \
	bcftools norm \
	--threads ${task.cpus} \
	--multiallelics -indels \
	--output-type z \
	- 2>"${varscan_indel_multiallelics_stats}" \
	| \
	bcftools norm \
	--threads ${task.cpus} \
	--fasta-ref "${reference_genome_fasta_forVarscanBcftoolsNorm}" \
	--output-type z \
	--output "${final_varscan_indel_vcf}" \
	- 2>"${varscan_indel_realign_normalize_stats}"

	tabix "${final_varscan_indel_vcf}"
	"""
}

// END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \\


// ~~~~~~~~~~~~~~~~ MuTect2 ~~~~~~~~~~~~~~~~ \\
// START

// BCFtools Concat ~ prepare the gnomAD allele frequency reference VCF for MuTect2 process, if needed
process mutect2GnomadReferenceVcfPrep_bcftools {
	publishDir "references/hg38", mode: 'copy'

	input:
	path gnomad_ref_vcf_chromosomes1_9
	path gnomad_ref_vcf_chromosomes1_9_index
	path gnomad_ref_vcf_chromosomes10_22
	path gnomad_ref_vcf_chromosomesXYM_alts
	path gnomad_ref_vcf_chromosomesXYM_alts_index

	output:
	tuple path(mutect_gnomad_ref_vcf), path(mutect_gnomad_ref_vcf_index) into mutect_gnomad_ref_vcf_fromProcess

	when:
	params.mutect == "on" && params.mutect_ref_vcf_concatenated == "no"

	script:
	mutect_gnomad_ref_vcf = "af-only-gnomad.hg38.vcf.gz"
	mutect_gnomad_ref_vcf_index = "${mutect_gnomad_ref_vcf}.tbi"
	"""
	bcftools concat \
	--threads ${task.cpus} \
	--output-type z \
	--output "${mutect_gnomad_ref_vcf}" \
	"${gnomad_ref_vcf_chromosomes1_9}" \
	"${gnomad_ref_vcf_chromosomes10_22}" \
	"${gnomad_ref_vcf_chromosomesXYM_alts}"

	tabix "${mutect_gnomad_ref_vcf}"
	"""
}

// Depending on whether the gnomAD allele frequency reference VCF was pre-built, set the input
// channel for the for MuTect2 process
if( params.mutect_ref_vcf_concatenated == "yes" && params.mutect == "on") {
	mutect_gnomad_ref_vcf = mutect_gnomad_ref_vcf_preBuilt.combine( mutect_gnomad_ref_vcf_index_preBuilt )
}
else {
	mutect_gnomad_ref_vcf = mutect_gnomad_ref_vcf_fromProcess
}

// Combine all reference FASTA files, WGS BED file, and resource VCFs into one channel for use in MuTect calling process
reference_genome_fasta_forMutectCalling.combine( reference_genome_fasta_index_forMutectCalling )
	.combine( reference_genome_fasta_dict_forMutectCalling )
	.combine( gatk_bundle_wgs_bed_forMutectCalling )
	.set{ reference_genome_bundle_and_bed_forMutectCalling }

reference_genome_bundle_and_bed_forMutectCalling.combine( mutect_gnomad_ref_vcf )
	.combine( panel_of_normals_1000G )
	.combine( panel_of_normals_1000G_index )
	.set{ reference_genome_bed_and_vcfs_forMutectCalling } 

// GATK MuTect2 ~ call somatic SNVs and indels via local assembly of haplotypes
process snvAndIndelCalling_gatk {
	publishDir "${params.output_dir}/somatic/mutect/rawVcfs", mode: 'symlink'
	tag "C=${chromosome} T=${tumor_id} N=${normal_id}"

	input:
	tuple path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(reference_genome_fasta_forMutectCalling), path(reference_genome_fasta_index_forMutectCalling), path(reference_genome_fasta_dict_forMutectCalling), path(gatk_bundle_wgs_bed_forMutectCalling), path(mutect_gnomad_ref_vcf), path(mutect_gnomad_ref_vcf_index), path(panel_of_normals_1000G), path(panel_of_normals_1000G_index) from tumor_normal_pair_forMutectCalling.combine(reference_genome_bed_and_vcfs_forMutectCalling)
	each chromosome from chromosome_list_forMutectCalling

	output:
	tuple val(tumor_normal_sample_id), path(raw_per_chromosome_vcf), path(raw_per_chromosome_vcf_index) into raw_per_chromosome_vcfs_forMutectVcfMerge
	tuple val(tumor_normal_sample_id), path(raw_per_chromosome_mutect_stats_file) into raw_per_chromosome_mutect_stats_forMutectStatsMerge

	when:
	params.mutect == "on"

	script:
	tumor_id = "${tumor_bam.baseName}".replaceFirst(/\..*$/, "")
	normal_id = "${normal_bam.baseName}".replaceFirst(/\..*$/, "")
	tumor_normal_sample_id = "${tumor_id}_vs_${normal_id}"
	per_chromosome_bed_file = "${gatk_bundle_wgs_bed_forMutectCalling}".replaceFirst(/\.bed/, ".${chromosome}.bed")
	raw_per_chromosome_vcf = "${tumor_normal_sample_id}.${chromosome}.vcf.gz"
	raw_per_chromosome_vcf_index = "${raw_per_chromosome_vcf}.tbi"
	raw_per_chromosome_mutect_stats_file = "${tumor_normal_sample_id}.${chromosome}.vcf.gz.stats"
	"""
	grep -w '${chromosome}' "${gatk_bundle_wgs_bed_forMutectCalling}" > "${per_chromosome_bed_file}"

	gatk Mutect2 \
	--java-options "-Xmx${task.memory.toGiga()}G -Djava.io.tmpdir=." \
	--verbosity ERROR \
	--tmp-dir . \
	--native-pair-hmm-threads ${task.cpus} \
	--af-of-alleles-not-in-resource 0.00003125 \
	--seconds-between-progress-updates 600 \
	--reference "${reference_genome_fasta_forMutectCalling}" \
	--intervals "${per_chromosome_bed_file}" \
	--germline-resource "${mutect_gnomad_ref_vcf}" \
	--panel-of-normals "${panel_of_normals_1000G}" \
	--input "${tumor_bam}" \
	--input "${normal_bam}" \
	--normal-sample "${normal_id}" \
	--output "${raw_per_chromosome_vcf}"
	"""
}

// GATK SortVcfs ~ merge all per chromosome MuTect2 VCFs
process mergeAndSortMutect2Vcfs_gatk {
	publishDir "${params.output_dir}/somatic/mutect/rawVcfs", mode: 'symlink'
	tag "${tumor_normal_sample_id}"
	
	input:
	tuple val(tumor_normal_sample_id), path(raw_per_chromosome_vcf), path(raw_per_chromosome_vcf_index) from raw_per_chromosome_vcfs_forMutectVcfMerge.groupTuple()

	output:
	tuple val(tumor_normal_sample_id), path(merged_raw_vcf), path(merged_raw_vcf_index) into merged_raw_vcfs_forMutectFilter

	when:
	params.mutect == "on"

	script:
	merged_raw_vcf = "${tumor_normal_sample_id}.vcf.gz"
	merged_raw_vcf_index = "${merged_raw_vcf}.tbi"
	"""
	gatk SortVcf \
	--java-options "-Xmx${task.memory.toGiga()}G -Djava.io.tmpdir=." \
	--VERBOSITY ERROR \
	--TMP_DIR . \
	--MAX_RECORDS_IN_RAM 4000000 \
	${raw_per_chromosome_vcf.collect { "--INPUT $it " }.join()} \
	--OUTPUT "${merged_raw_vcf}"
	"""
}

// GATK MergeMutectStats ~ merge the per chromosome stats output of MuTect2 variant calling
process mergeMutect2StatsForFiltering_gatk {
	publishDir "${params.output_dir}/somatic/mutect/rawVcfs", mode: 'symlink'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(raw_per_chromosome_mutect_stats_file) from raw_per_chromosome_mutect_stats_forMutectStatsMerge.groupTuple()

	output:
	tuple val(tumor_normal_sample_id), path(merged_mutect_stats_file) into merged_mutect_stats_file_forMutectFilter

	when:
	params.mutect == "on"

	script:
	merged_mutect_stats_file = "${tumor_normal_sample_id}.vcf.gz.stats"
	"""
	gatk MergeMutectStats \
	--java-options "-Xmx${task.memory.toGiga()}G -Djava.io.tmpdir=." \
	--verbosity ERROR \
	--tmp-dir . \
	${raw_per_chromosome_mutect_stats_file.collect { "--stats $it " }.join()} \
	--output "${merged_mutect_stats_file}"
	"""
}

// Combine WGS BED file and resource VCFs into one channel for use in MuTect pileup process
gatk_bundle_wgs_bed_forMutectPileup.combine( exac_common_sites_ref_vcf )
	.combine( exac_common_sites_ref_vcf_index )
	.set{ bed_and_resources_vcfs_forMutectPileup }

// GATK GetPileupSummaries ~ tabulates pileup metrics for inferring contamination
process pileupSummariesForMutect2Contamination_gatk {
	publishDir "${params.output_dir}/somatic/mutect/pileups", mode: 'symlink'
	tag "C=${chromosome} T=${tumor_id} N=${normal_id}"

	input:
	tuple path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(gatk_bundle_wgs_bed_forMutectPileup), path(exac_common_sites_ref_vcf), path(exac_common_sites_ref_vcf_index) from tumor_normal_pair_forMutectPileup.combine(bed_and_resources_vcfs_forMutectPileup)
	each chromosome from chromosome_list_forMutectPileup

	output:
	tuple val(tumor_normal_sample_id), path(per_chromosome_tumor_pileup) into per_chromosome_tumor_pileups_forMutectPileupGather
	tuple val(tumor_normal_sample_id), path(per_chromosome_normal_pileup) into per_chromosome_normal_pileups_forMutectPileupGather

	when:
	params.mutect == "on"

	script:
	tumor_id = "${tumor_bam.baseName}".replaceFirst(/\..*$/, "")
	normal_id = "${normal_bam.baseName}".replaceFirst(/\..*$/, "")
	tumor_normal_sample_id = "${tumor_id}_vs_${normal_id}"
	per_chromosome_bed_file = "${gatk_bundle_wgs_bed_forMutectPileup}".replaceFirst(/\.bed/, ".${chromosome}.bed")
	per_chromosome_tumor_pileup = "${tumor_id}.${chromosome}.pileup"
	per_chromosome_normal_pileup = "${normal_id}.${chromosome}.pileup"
	"""
	grep -w '${chromosome}' "${gatk_bundle_wgs_bed_forMutectPileup}" > "${per_chromosome_bed_file}"

	gatk GetPileupSummaries \
	--java-options "-Xmx${task.memory.toGiga()}G -Djava.io.tmpdir=." \
	--verbosity ERROR \
	--tmp-dir . \
	--intervals "${per_chromosome_bed_file}" \
	--variant "${exac_common_sites_ref_vcf}" \
	--input "${tumor_bam}" \
	--output "${per_chromosome_tumor_pileup}"

	gatk GetPileupSummaries \
	--java-options "-Xmx${task.memory.toGiga()}G -Djava.io.tmpdir=." \
	--verbosity ERROR \
	--tmp-dir . \
	--intervals "${per_chromosome_bed_file}" \
	--variant "${exac_common_sites_ref_vcf}" \
	--input "${normal_bam}" \
	--output "${per_chromosome_normal_pileup}"
	"""
}

// GATK GatherPileupSummaries ~ combine tumor pileup tables for inferring contamination
process gatherTumorPileupSummariesForMutect2Contamination_gatk {
	publishDir "${params.output_dir}/somatic/mutect/pileups", mode: 'symlink'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(per_chromosome_tumor_pileup), path(reference_genome_fasta_dict) from per_chromosome_tumor_pileups_forMutectPileupGather.groupTuple().combine(reference_genome_fasta_dict_forMutectPileupGatherTumor)

	output:
	tuple val(tumor_normal_sample_id), path(tumor_pileup) into tumor_pileups_forMutectContamination

	when:
	params.mutect == "on"

	script:
	tumor_id = "${tumor_normal_sample_id}".replaceFirst(/\_vs\_.*/, "")
	tumor_pileup = "${tumor_id}.pileup"
	"""
	gatk GatherPileupSummaries \
	--java-options "-Xmx${task.memory.toGiga()}G -Djava.io.tmpdir=." \
	--verbosity ERROR \
	--tmp-dir . \
	--sequence-dictionary "${reference_genome_fasta_dict}" \
	${per_chromosome_tumor_pileup.collect { "--I $it " }.join()} \
	--O "${tumor_pileup}"
	"""
}

// GATK GatherPileupSummaries ~ combine normal pileup tables for inferring contamination
process gatherNormalPileupSummariesForMutect2Contamination_gatk {
	publishDir "${params.output_dir}/somatic/mutect/pileups", mode: 'symlink'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(per_chromosome_normal_pileup), path(reference_genome_fasta_dict) from per_chromosome_normal_pileups_forMutectPileupGather.groupTuple().combine(reference_genome_fasta_dict_forMutectPileupGatherNormal)

	output:
	tuple val(tumor_normal_sample_id), path(normal_pileup) into normal_pileups_forMutectContamination

	when:
	params.mutect == "on"

	script:
	normal_id = "${tumor_normal_sample_id}".replaceFirst(/.*\_vs\_/, "")
	normal_pileup = "${normal_id}.pileup"
	"""
	gatk GatherPileupSummaries \
	--java-options "-Xmx${task.memory.toGiga()}G -Djava.io.tmpdir=." \
	--verbosity ERROR \
	--tmp-dir . \
	--sequence-dictionary "${reference_genome_fasta_dict}" \
	${per_chromosome_normal_pileup.collect { "--I $it " }.join()} \
	--O "${normal_pileup}"
	"""
}

// GATK CalculateContamination ~ calculate the fraction of reads coming from cross-sample contamination
process mutect2ContaminationCalculation_gatk {
	publishDir "${params.output_dir}/somatic/mutect", mode: 'copy'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(tumor_pileup), path(normal_pileup) from tumor_pileups_forMutectContamination.join(normal_pileups_forMutectContamination)

	output:
	tuple val(tumor_normal_sample_id), path(contamination_file) into contamination_file_forMutectFilter

	when:
	params.mutect == "on"

	script:
	contamination_file = "${tumor_normal_sample_id}.mutect.contamination.txt" 
	"""
	gatk CalculateContamination \
	--java-options "-Xmx${task.memory.toGiga()}G -Djava.io.tmpdir=." \
	--verbosity ERROR \
	--tmp-dir . \
	--input "${tumor_pileup}" \
	--matched-normal "${normal_pileup}" \
	--output "${contamination_file}"
	"""
}

// Combine all reference FASTA files and input VCF, stats file, and contamination table into one channel for use in MuTect filtering process
reference_genome_fasta_forMutectFilter.combine( reference_genome_fasta_index_forMutectFilter )
	.combine( reference_genome_fasta_dict_forMutectFilter )
	.set{ reference_genome_bundle_forMutectFilter }

merged_raw_vcfs_forMutectFilter.join( merged_mutect_stats_file_forMutectFilter )
	.join( contamination_file_forMutectFilter )
	.set{ input_vcf_stats_and_contamination_forMutectFilter }

// GATK FilterMutectCalls ~ filter somatic SNVs and indels called by Mutect2
process mutect2VariantFiltration_gatk {
	publishDir "${params.output_dir}/somatic/mutect/filteredVcfs", mode: 'symlink'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(merged_raw_vcf), path(merged_raw_vcf_index), path(merged_mutect_stats_file), path(contamination_file), path(reference_genome_fasta_forMutectFilter), path(reference_genome_fasta_index_forMutectFilter), path(reference_genome_fasta_dict_forMutectFilter) from input_vcf_stats_and_contamination_forMutectFilter.combine(reference_genome_bundle_forMutectFilter)

	output:
	tuple val(tumor_normal_sample_id), path(filtered_vcf), path(filtered_vcf_index) into filtered_vcf_forMutectBcftools
	path filter_stats_file

	when:
	params.mutect == "on"

	script:
	filtered_vcf = "${tumor_normal_sample_id}.filtered.vcf.gz"
	filtered_vcf_index = "${filtered_vcf}.tbi"
	filter_stats_file = "${tumor_normal_sample_id}.filterstats.txt"
	"""
	gatk FilterMutectCalls \
	--java-options "-Xmx${task.memory.toGiga()}G -Djava.io.tmpdir=." \
	--verbosity ERROR \
	--tmp-dir . \
	--unique-alt-read-count 5 \
	--reference "${reference_genome_fasta_forMutectFilter}" \
	--stats "${merged_mutect_stats_file}" \
	--variant "${merged_raw_vcf}" \
	--contamination-table "${contamination_file}" \
	--output "${filtered_vcf}" \
	--filtering-stats "${filter_stats_file}"
	"""
}

// Combine all needed reference FASTA files into one channel for use in MuTect2 / BCFtools Norm process
reference_genome_fasta_forMutectBcftools.combine( reference_genome_fasta_index_forMutectBcftools )
	.combine( reference_genome_fasta_dict_forMutectBcftools )
	.set{ reference_genome_bundle_forMutectBcftools }

// BCFtools Norm ~ split multiallelic sites into multiple rows then left-align and normalize indels
process splitMultiallelicAndLeftNormalizeMutect2Vcf_bcftools {
	publishDir "${params.output_dir}/somatic/mutect", mode: 'copy'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(filtered_vcf), path(filtered_vcf_index), path(reference_genome_fasta_forMutectBcftools), path(reference_genome_fasta_index_forMutectBcftools), path(reference_genome_fasta_dict_forMutectBcftools) from filtered_vcf_forMutectBcftools.combine(reference_genome_bundle_forMutectBcftools)

	output:
	path final_mutect_vcf into final_mutect_vcf_forAnnotation
	path final_mutect_vcf_index into final_mutect_vcf_index_forAnnotation
	tuple val(tumor_normal_sample_id), path(final_mutect_vcf), path(final_mutect_vcf_index) into mutect_vcf_forSclust
	tuple val(tumor_normal_sample_id), path(final_mutect_vcf), path(final_mutect_vcf_index) into final_combined_mutect_vcf
	path mutect_multiallelics_stats
	path mutect_realign_normalize_stats

	when:
	params.mutect == "on"

	script:
	final_mutect_vcf = "${tumor_normal_sample_id}.mutect.somatic.vcf.gz"
	final_mutect_vcf_index = "${final_mutect_vcf}.tbi"
	mutect_multiallelics_stats = "${tumor_normal_sample_id}.mutect.multiallelicsstats.txt"
	mutect_realign_normalize_stats = "${tumor_normal_sample_id}.mutect.realignnormalizestats.txt"
	"""
	zcat "${filtered_vcf}" \
	| \
	grep -E '^#|PASS' \
	| \
	bcftools norm \
	--threads ${task.cpus} \
	--multiallelics -both \
	--output-type z \
	- 2>"${mutect_multiallelics_stats}" \
	| \
	bcftools norm \
	--threads ${task.cpus} \
	--fasta-ref "${reference_genome_fasta_forMutectBcftools}" \
	--output-type z \
	- 2>"${mutect_realign_normalize_stats}" \
	--output "${final_mutect_vcf}"

	tabix "${final_mutect_vcf}"
	"""
}

// BCFtools view ~ separate out normalized SNVs and indel calls for consensus call generation
process splitMutectSnvsAndIndelsForConsensus_bcftools {
	publishDir "${params.output_dir}/somatic/mutect", mode: 'copy', pattern: '*.{vcf.gz,tbi}'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(final_mutect_vcf), path(final_mutect_vcf_index) from final_combined_mutect_vcf

	output:
	tuple val(tumor_normal_sample_id), path(final_mutect_snv_vcf), path(final_mutect_snv_vcf_index) into final_mutect_snv_vcf_forConsensus
	tuple val(tumor_normal_sample_id), path(final_mutect_indel_vcf), path(final_mutect_indel_vcf_index) into final_mutect_indel_vcf_forConsensus

	when:
	params.mutect == "on"

	script:
	final_mutect_snv_vcf = "${final_mutect_vcf}".replaceFirst(/\.vcf\.gz/, ".snv.vcf.gz")
	final_mutect_snv_vcf_index ="${final_mutect_snv_vcf}.tbi"
	final_mutect_indel_vcf = "${final_mutect_vcf}".replaceFirst(/\.vcf\.gz/, ".indel.vcf.gz")
	final_mutect_indel_vcf_index = "${final_mutect_indel_vcf}.tbi"
	"""
	bcftools view \
	--threads "${task.cpus}" \
	--output-type z \
	--types snps \
	--output-file "${final_mutect_snv_vcf}" \
	"${final_mutect_vcf}"

	tabix "${final_mutect_snv_vcf}"

	bcftools view \
	--threads "${task.cpus}" \
	--output-type z \
	--types indels \
	--output-file "${final_mutect_indel_vcf}" \
	"${final_mutect_vcf}"

	tabix "${final_mutect_indel_vcf}"
	"""
}

// END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \\


// ~~~~~~~~~~~~~~~~ ascatNGS ~~~~~~~~~~~~~~~ \\
// START

// Combine all reference files into one channel for us in AscatNGS
reference_genome_fasta_forAscatNgs.combine( reference_genome_fasta_index_forAscatNgs )
	.combine( reference_genome_fasta_dict_forAscatNgs )
	.combine( snp_gc_corrections )
	.set{ reference_genome_and_snpgc_forAscatNgs }

// ascatNGS ~  identifying somatically acquired copy-number alterations
process cnvCalling_ascatngs {
	publishDir "${params.output_dir}/somatic/ascatNgs", mode: 'copy', pattern: '*.{png,csv,vcf.gz,txt,tbi}'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(sample_sex), path(reference_genome_fasta_forAscatNgs), path(reference_genome_fasta_index_forAscatNgs), path(reference_genome_fasta_dict_forAscatNgs), path(snp_gc_corrections) from bams_and_sex_of_sample_forAscatNgs.combine(reference_genome_and_snpgc_forAscatNgs)

	output:
	tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(cnv_profile_final), path(run_statistics) into bams_cnv_profile_and_statistics_forCaveman
	tuple path(cnv_profile_vcf), path(cnv_profile_vcf_index)
	path ascat_profile_png
	path ascat_raw_profile_png
	path sunrise_png
	path aspcf_png
	path germline_png
	path tumor_png

	when:
	params.ascatngs == "on"

	script:
	ploidy_and_purity = (params.ascatngs_ploidy && params.ascatngs_purity) ? "-ploidy ${params.ascatngs_ploidy} -purity ${params.ascatngs_purity}" : ""
	tumor_id = "${tumor_bam.baseName}".replaceFirst(/\..*$/, "")
	ascat_profile_png = "${tumor_normal_sample_id}.ascat.profile.png"
	ascat_raw_profile_png = "${tumor_normal_sample_id}.ascat.raw.profile.png"
	sunrise_png = "${tumor_normal_sample_id}.ascat.sunrise.png"
	aspcf_png = "${tumor_normal_sample_id}.ascat.aspcf.png"
	germline_png = "${tumor_normal_sample_id}.ascat.germline.png"
	tumor_png = "${tumor_normal_sample_id}.ascat.tumor.png"
	cnv_profile_final = "${tumor_normal_sample_id}.ascat.cnv.csv"
	cnv_profile_vcf = "${tumor_normal_sample_id}.ascat.vcf.gz"
	cnv_profile_vcf_index = "${cnv_profile_vcf}.tbi"
	run_statistics = "${tumor_normal_sample_id}.ascat.runstatistics.txt"
	"""
	sex=\$(cut -d ' ' -f 2 "${sample_sex}")

	ascat.pl \
	-outdir . \
	-tumour "${tumor_bam}" \
	-normal "${normal_bam}" \
	-reference "${reference_genome_fasta_forAscatNgs}" \
	-snp_gc "${snp_gc_corrections}" \
	-protocol WGS \
	-gender "\${sex}" \
	-genderChr chrY \
	-species Homo_sapiens \
	-assembly GRCh38 \
	-cpus "${task.cpus}" \
	-nobigwig \
	${ploidy_and_purity}

	mv "${tumor_id}.ASCATprofile.png" "${ascat_profile_png}"
	mv "${tumor_id}.rawprofile.png" "${ascat_raw_profile_png}"
	mv "${tumor_id}.sunrise.png" "${sunrise_png}"
	mv "${tumor_id}.ASPCF.png" "${aspcf_png}"
	mv "${tumor_id}.germline.png" "${germline_png}"
	mv "${tumor_id}.tumour.png" "${tumor_png}"
	mv "${tumor_id}.copynumber.caveman.vcf.gz" "${cnv_profile_vcf}"
	mv "${tumor_id}.copynumber.caveman.vcf.gz.tbi" "${cnv_profile_vcf_index}"
	mv "${tumor_id}.samplestatistics.txt" "${run_statistics}"

	touch "${cnv_profile_final}"
	echo "segment_number,chromosome,start_position,end_position,normal_total_copy_number,normal_minor_copy_number,tumor_total_copy_number,tumor_minor_copy_number" >> "${cnv_profile_final}"
	cat "${tumor_id}.copynumber.caveman.csv" >> "${cnv_profile_final}"
	"""
}

// END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \\


// ~~~~~~~~~~~~~ Control-FREEC ~~~~~~~~~~~~~ \\
// START

// Combine all reference FASTA files and WGS BED file into one channel for use in Control-FREEC / SAMtools mpileup
reference_genome_fasta_forControlFreecSamtoolsMpileup.combine( reference_genome_fasta_index_forControlFreecSamtoolsMpileup )
	.combine( reference_genome_fasta_dict_forControlFreecSamtoolsMpileup )
	.set{ reference_genome_bundle_forControlFreecSamtoolsMpileup }

reference_genome_bundle_forControlFreecSamtoolsMpileup.combine( gatk_bundle_wgs_bed_forControlFreecSamtoolsMpileup )
	.set{ reference_genome_bundle_and_bed_forControlFreecSamtoolsMpileup }

// SAMtools mpileup ~ generate per chromosome mpileups for the tumor and normal BAMs separately
process bamMpileupForControlFreec_samtools {
	publishDir "${params.output_dir}/somatic/controlFreec/pileups", mode: 'symlink'
	tag "C=${chromosome} T=${tumor_id} N=${normal_id}"

	input:
	tuple path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(reference_genome_fasta_forControlFreecSamtoolsMpileup), path(reference_genome_fasta_index_forControlFreecSamtoolsMpileup), path(reference_genome_fasta_dict_forControlFreecSamtoolsMpileup), path(gatk_bundle_wgs_bed_forControlFreecSamtoolsMpileup) from tumor_normal_pair_forControlFreecSamtoolsMpileup.combine(reference_genome_bundle_and_bed_forControlFreecSamtoolsMpileup)
	each chromosome from chromosome_list_forControlFreecSamtoolsMpileup

	output:
	tuple val(tumor_normal_sample_id), val(tumor_id), val(normal_id), path(tumor_pileup_per_chromosome), path(normal_pileup_per_chromosome) into per_chromosome_tumor_normal_pileups_forControlFreecMerge

	when:
	params.controlfreec == "on"

	script:
	tumor_id = "${tumor_bam.baseName}".replaceFirst(/\..*$/, "")
	normal_id = "${normal_bam.baseName}".replaceFirst(/\..*$/, "")
	tumor_normal_sample_id = "${tumor_id}_vs_${normal_id}"
	tumor_pileup_per_chromosome = "${tumor_id}.${chromosome}.pileup.gz"
	normal_pileup_per_chromosome = "${normal_id}.${chromosome}.pileup.gz"
	"""
	samtools mpileup \
	--positions "${gatk_bundle_wgs_bed_forControlFreecSamtoolsMpileup}" \
	--fasta-ref "${reference_genome_fasta_forControlFreecSamtoolsMpileup}" \
	--region "${chromosome}" \
	"${tumor_bam}" \
	| \
	bgzip > "${tumor_pileup_per_chromosome}"

	samtools mpileup \
	--positions "${gatk_bundle_wgs_bed_forControlFreecSamtoolsMpileup}" \
	--fasta-ref "${reference_genome_fasta_forControlFreecSamtoolsMpileup}" \
	--region "${chromosome}" \
	"${normal_bam}" \
	| \
	bgzip > "${normal_pileup_per_chromosome}"
	"""
}

// Mpileup Merge ~ Combine all tumor and normal per chromosome mpileup files separately
process mergeMpileupsForControlFreec_samtools {
	publishDir "${params.output_dir}/somatic/controlFreec/pileups", mode: 'symlink'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), val(tumor_id), val(normal_id), path(tumor_pileup_per_chromosome), path(normal_pileup_per_chromosome) from per_chromosome_tumor_normal_pileups_forControlFreecMerge.groupTuple(by: [0,1,2])
	val chromosome_list from chromosome_list_forControlFreecMerge.collect()

	output:
	tuple val(tumor_normal_sample_id), path(tumor_pileup), path(normal_pileup) into tumor_normal_pileups_forControlFreecCalling

	when:
	params.controlfreec == "on"

	script:
	tumor_pileup = "${tumor_id}.pileup.gz"
	normal_pileup = "${normal_id}.pileup.gz"
	"""
	chromosome_array=\$(echo "${chromosome_list}" | sed 's/[][]//g' | sed 's/,//g')

	for chrom in \$chromosome_array;
		do
			zcat "${tumor_id}.\${chrom}.pileup.gz" >> "${tumor_id}.pileup"

			zcat "${normal_id}.\${chrom}.pileup.gz" >> "${normal_id}.pileup"
		done

	bgzip < "${tumor_id}.pileup" > "${tumor_pileup}"
	bgzip < "${normal_id}.pileup" > "${normal_pileup}"
	"""
}

// Combine mpileup input files with the sample sex identity then all reference files into one channel for use in Control-FREEC
tumor_normal_pileups_forControlFreecCalling.join( sex_of_sample_forControlFreecCalling )
	.set{ tumor_normal_pileups_and_sex_ident }

reference_genome_fasta_forControlFreecCalling.combine( reference_genome_fasta_index_forControlFreecCalling )
	.combine( reference_genome_fasta_dict_forControlFreecCalling )
	.set{ reference_genome_bundle_forControlFreecCalling }

reference_genome_bundle_forControlFreecCalling.combine( autosome_sex_chromosome_fasta_dir )
	.combine( autosome_sex_chromosome_sizes )
	.combine( common_dbsnp_ref_vcf )
	.combine( common_dbsnp_ref_vcf_index )
	.combine( mappability_track_zip )
	.set{ reference_data_bundle_forControlFreec }

// Control-FREEC ~ detection of copy-number changes and allelic imbalances
process cnvCalling_controlfreec {
	publishDir "${params.output_dir}/somatic/controlFreec", mode: 'copy', pattern: '*.{txt,cnv}'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(tumor_pileup), path(normal_pileup), path(sex_of_sample_forControlFreecCalling), path(reference_genome_fasta_forControlFreecCalling), path(reference_genome_fasta_index_forControlFreecCalling), path(reference_genome_fasta_dict_forControlFreecCalling), path(autosome_sex_chromosome_fasta_dir), path(autosome_sex_chromosome_sizes), path(common_dbsnp_ref_vcf), path(common_dbsnp_ref_vcf_index), path(mappability_track_zip) from tumor_normal_pileups_and_sex_ident.combine(reference_data_bundle_forControlFreec)

	output:
	tuple val(tumor_normal_sample_id), path(cnv_profile_raw), path(cnv_ratio_file), path(baf_file) into cnv_calling_files_forControlFreecPostProcessing
	path control_freec_config_file
	path subclones_file

	when:
	params.controlfreec == "on"

	script:
	control_freec_config_file = "${tumor_normal_sample_id}.controlfreec.config.txt"
	cnv_profile_raw = "${tumor_normal_sample_id}.controlfreec.raw.cnv"
	cnv_ratio_file = "${tumor_normal_sample_id}.controlfreec.ratio.txt"
	subclones_file = "${tumor_normal_sample_id}.controlfreec.subclones.txt"
	baf_file = "${tumor_normal_sample_id}.controlfreec.baf.txt"
	"""
	unzip -q "${mappability_track_zip}"
	sex=\$(cut -d ' ' -f 2 "${sex_of_sample_forControlFreecCalling}")

	touch "${control_freec_config_file}"
	echo "[general]" >> "${control_freec_config_file}"
	echo "chrFiles = \${PWD}/${autosome_sex_chromosome_fasta_dir}" >> "${control_freec_config_file}"
	echo "chrLenFile = \${PWD}/${autosome_sex_chromosome_sizes}" >> "${control_freec_config_file}"
	echo "gemMappabilityFile = \${PWD}/out100m2_hg38.gem" >> "${control_freec_config_file}"
	echo "minimalSubclonePresence = 20" >> "${control_freec_config_file}"
	echo "maxThreads = ${task.cpus}" >> "${control_freec_config_file}"
	echo "ploidy = 2" "${control_freec_config_file}"
	echo "sex = \${sex}" >> "${control_freec_config_file}"
	echo "window = 50000" >> "${control_freec_config_file}"
	echo "" >> "${control_freec_config_file}"

	echo "[sample]" >> "${control_freec_config_file}"
	echo "mateFile = ${tumor_pileup}" >> "${control_freec_config_file}"
	echo "inputFormat = pileup" >> "${control_freec_config_file}"
	echo "mateOrientation = FR" >> "${control_freec_config_file}"
	echo "" >> "${control_freec_config_file}"

	echo "[control]" >> "${control_freec_config_file}"
	echo "mateFile = ${normal_pileup}" >> "${control_freec_config_file}"
	echo "inputFormat = pileup" >> "${control_freec_config_file}"
	echo "mateOrientation = FR" >> "${control_freec_config_file}"
	echo "" >> "${control_freec_config_file}"

	echo "[BAF]" >> "${control_freec_config_file}"
	echo "SNPfile = \${PWD}/${common_dbsnp_ref_vcf}" >> "${control_freec_config_file}"
	echo "" >> "${control_freec_config_file}"

	freec -conf "${control_freec_config_file}"

	mv "${tumor_pileup}_CNVs" "${cnv_profile_raw}"
	mv "${tumor_pileup}_ratio.txt" "${cnv_ratio_file}"
	mv "${tumor_pileup}_subclones.txt" "${subclones_file}"
	mv "${tumor_pileup}_BAF.txt" "${baf_file}"
	"""
}

// Control-FREEC ~ post-processing of CNV predictions for significance, visualization, and format compatibility
process cnvPredictionPostProcessing_controlfreec {
	publishDir "${params.output_dir}/somatic/controlFreec", mode: 'copy', pattern: '*.{txt,bed,png}'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(cnv_profile_raw), path(cnv_ratio_file), path(baf_file) from cnv_calling_files_forControlFreecPostProcessing

	output:
	path cnv_profile_final
	path cnv_ratio_bed_file
	path ratio_graph_png
	path ratio_log2_graph_png
	path baf_graph_png

	when:
	params.controlfreec == "on"

	script:
	cnv_profile_final = "${tumor_normal_sample_id}.controlfreec.cnv.txt"
	cnv_ratio_bed_file = "${tumor_normal_sample_id}.controlfreec.ratio.bed"
	ratio_graph_png = "${tumor_normal_sample_id}.controlfreec.ratio.png"
	ratio_log2_graph_png = "${tumor_normal_sample_id}.controlfreec.ratio.log2.png"
	baf_graph_png = "${tumor_normal_sample_id}.controlfreec.baf.png"
	"""
	cat \${CONTROLFREEC_DIR}/scripts/assess_significance.R | R --slave --args "${cnv_profile_raw}" "${cnv_ratio_file}"
	mv "${cnv_profile_raw}.p.value.txt" "${cnv_profile_final}"

	cat \${CONTROLFREEC_DIR}/scripts/makeGraph.R | R --slave --args 2 "${cnv_ratio_file}" "${baf_file}"
	mv "${cnv_ratio_file}.png" "${ratio_graph_png}"
	mv "${cnv_ratio_file}.log2.png" "${ratio_log2_graph_png}"
	mv "${baf_file}.png" "${baf_graph_png}"

	perl \${CONTROLFREEC_DIR}/scripts/freec2bed.pl -f "${cnv_ratio_file}" > "${cnv_ratio_bed_file}"
	"""
}

// END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \\


// ~~~~~~~~~~~~~~~~~ Sclust ~~~~~~~~~~~~~~~~~ \\
// START

// Sclust bamprocess ~ extract the read ratio and SNP information per chromosome
process bamprocessPerChromosome_sclust {
	publishDir "${params.output_dir}/somatic/sclust/intermediates", mode: 'symlink'
	tag "C=${chromosome} T=${tumor_id} N=${normal_id}"

	input:
	tuple path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index) from tumor_normal_pair_forSclustBamprocess
	each chromosome from chromosome_list_forSclustBamprocess

	output:
	tuple val(tumor_normal_sample_id), path(bamprocess_data_per_chromosome) into per_chromosome_bamprocess_data

	when:
	params.sclust == "on" && params.mutect == "on"

	script:
	tumor_id = "${tumor_bam.baseName}".replaceFirst(/\..*$/, "")
	normal_id = "${normal_bam.baseName}".replaceFirst(/\..*$/, "")
	tumor_normal_sample_id = "${tumor_id}_vs_${normal_id}"
	bamprocess_data_per_chromosome = "${tumor_normal_sample_id}_${chromosome}_bamprocess_data.txt"
	"""
	Sclust bamprocess \
	-t "${tumor_bam}" \
	-n "${normal_bam}" \
	-o "${tumor_normal_sample_id}" \
	-build hg38 \
	-part 2 \
	-r "${chromosome}"
	"""
}

// Sclust bamprocess ~ merge per chromosome bamprocess data files and generate a read-count and common SNP-count files
process mergeBamprocessData_sclust {
	publishDir "${params.output_dir}/somatic/sclust/intermediates", mode: 'symlink', pattern: '*.{txt}'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(bamprocess_data_per_chromosome) from per_chromosome_bamprocess_data.groupTuple()

	output:
	tuple val(tumor_normal_sample_id), path(read_count_file), path(common_snp_count_file) into read_count_and_snp_count_files

	when:
	params.sclust == "on" && params.mutect == "on"

	script:
	read_count_file = "${tumor_normal_sample_id}_rcount.txt"
	common_snp_count_file = "${tumor_normal_sample_id}_snps.txt"
	"""
	Sclust bamprocess \
	-build hg38 \
	-i "${tumor_normal_sample_id}" \
	-o "${tumor_normal_sample_id}"
	"""
}

// VCFtools / VAtools ~ extract metrics from VCF to prepare for use in Sclust process
process prepareVcfForSclust_vcftools {
	publishDir "${params.output_dir}/somatic/sclust/intermediates", mode: 'symlink', pattern: '*.{vcf.gz}'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(final_mutect_vcf), path(final_mutect_vcf_index) from mutect_vcf_forSclust

	output:
	tuple val(tumor_normal_sample_id), path(mutations_vcf) into vcf_forSclustCn

	when:
	params.sclust == "on" && params.mutect == "on"

	script:
	tumor_id = "${tumor_normal_sample_id}".replaceFirst(/\_vs\_.*/, "")
	normal_id = "${tumor_normal_sample_id}".replaceFirst(/.*\_vs\_/, "")
	tumor_normal_read_depth = "${tumor_normal_sample_id}.DP.FORMAT"
	tumor_normal_allelic_frequency = "${tumor_normal_sample_id}.AF.FORMAT"
	tumor_normal_forward_reverse_reads = "${tumor_normal_sample_id}.SB.FORMAT"
	mutations_vcf = "${tumor_normal_sample_id}.sclust.final.vcf.gz"
	"""
	zcat "${final_mutect_vcf}" \
	| \
	sed 's|Approximate read depth; some reads may have been filtered|Read Depth Tumor|' \
	| \
	gzip > "${tumor_normal_sample_id}.sclust.base.vcf.gz"

	vcftools \
	--gzvcf "${tumor_normal_sample_id}.sclust.base.vcf.gz" \
	--stdout \
	--extract-FORMAT-info DP \
	--indv "${tumor_id}" \
	| \
	grep -v 'CHROM' > tumor.dp.tsv

	vcftools \
	--gzvcf "${tumor_normal_sample_id}.sclust.base.vcf.gz" \
	--stdout \
	--extract-FORMAT-info DP \
	--indv "${normal_id}" \
	| \
	grep -v 'CHROM' > normal.dp.tsv

	vcftools \
	--gzvcf "${tumor_normal_sample_id}.sclust.base.vcf.gz" \
	--stdout \
	--extract-FORMAT-info AF \
	--indv "${tumor_id}" \
	| \
	grep -v 'CHROM' > tumor.af.tsv

	vcftools \
	--gzvcf "${tumor_normal_sample_id}.sclust.base.vcf.gz" \
	--stdout \
	--extract-FORMAT-info AF \
	--indv "${normal_id}" \
	| \
	grep -v 'CHROM' > normal.af.tsv

	vcftools \
	--gzvcf "${tumor_normal_sample_id}.sclust.base.vcf.gz" \
	--stdout \
	--extract-FORMAT-info SB \
	--indv "${tumor_id}" \
	| \
	grep -v 'CHROM' \
	| \
	forward_reverse_score_calculator.py > tumor.fr.tsv
	
	awk 'BEGIN {OFS="\t"} {print \$1,\$2,"."}' tumor.dp.tsv > placeholder.tg.tsv

	vcf-info-annotator \
	--overwrite \
	--description "Read Depth Tumor" \
	--value_format Integer \
	--output-vcf "${tumor_normal_sample_id}.sclust.int1.vcf.gz" \
	"${tumor_normal_sample_id}.sclust.base.vcf.gz" \
	tumor.dp.tsv \
	DP

	vcf-info-annotator \
	--description "Read Depth Normal" \
	--value_format Integer \
	--output-vcf "${tumor_normal_sample_id}.sclust.int2.vcf.gz" \
	"${tumor_normal_sample_id}.sclust.int1.vcf.gz" \
	normal.dp.tsv \
	DP_N

	vcf-info-annotator \
	--description "Allelic Frequency Tumor" \
	--value_format Float \
	--output-vcf "${tumor_normal_sample_id}.sclust.int3.vcf.gz" \
	"${tumor_normal_sample_id}.sclust.int2.vcf.gz" \
	tumor.af.tsv \
	AF

	vcf-info-annotator \
	--description "Allelic Frequency Normal" \
	--value_format Float \
	--output-vcf "${tumor_normal_sample_id}.sclust.int4.vcf.gz" \
	"${tumor_normal_sample_id}.sclust.int3.vcf.gz" \
	normal.af.tsv \
	AF_N

	vcf-info-annotator \
	--description "Forward-Reverse Score" \
	--value_format Float \
	--output-vcf "${tumor_normal_sample_id}.sclust.int5.vcf.gz" \
	"${tumor_normal_sample_id}.sclust.int4.vcf.gz" \
	tumor.fr.tsv \
	FR

	vcf-info-annotator \
	--description "Target Name (Genome Partition)" \
	--value_format String \
	--output-vcf "${mutations_vcf}" \
	"${tumor_normal_sample_id}.sclust.int5.vcf.gz" \
	placeholder.tg.tsv \
	TG
	"""
}

// Sclust cn / cluster ~ perform copy-number analysis, estimation of tumor purity, and mutational clustering
process cnvCalling_sclust {
	publishDir "${params.output_dir}/somatic/sclust", mode: 'copy', pattern: '*.{txt,pdf}'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(read_count_file), path(common_snp_count_file), path(mutations_vcf) from read_count_and_snp_count_files.join(vcf_forSclustCn)

	output:
	path mutations_exp_af_file
	path allelic_states_file
	path cnv_profile_pdf
	path cnv_profile_file
	path cnv_segments_file
	path sclust_subclones_file
	path mutation_clusters_file
	path mutation_clusters_pdf
	path cluster_assignment_file

	when:
	params.sclust == "on" && params.mutect == "on"

	script:
	allelic_states_file = "${tumor_normal_sample_id}.sclust.allelicstates.txt"
	cnv_profile_pdf = "${tumor_normal_sample_id}.sclust.profile.pdf"
	cnv_profile_file = "${tumor_normal_sample_id}.sclust.cnvsummary.txt"
	cnv_segments_file = "${tumor_normal_sample_id}.sclust.cnvsegments.txt"
	mutations_exp_af_file = "${tumor_normal_sample_id}_muts_expAF.txt"
	sclust_subclones_file = "${tumor_normal_sample_id}.sclust.subclones.txt"
	mutation_clusters_file = "${tumor_normal_sample_id}.sclust.mutclusters.txt"
	mutation_clusters_pdf = "${tumor_normal_sample_id}.sclust.mutclusters.pdf"
	cluster_assignment_file = "${tumor_normal_sample_id}.sclust.clusterassignments.txt"
	"""
	gunzip -f "${mutations_vcf}"

	Sclust cn \
	-ns 1000 \
	-rc "${read_count_file}" \
	-snp "${common_snp_count_file}" \
	-vcf "${tumor_normal_sample_id}.sclust.final.vcf" \
	-o "${tumor_normal_sample_id}"

	mv "${tumor_normal_sample_id}_allelic_states.txt" "${allelic_states_file}"
	mv "${tumor_normal_sample_id}_cn_profile.pdf" "${cnv_profile_pdf}"
	mv "${tumor_normal_sample_id}_cn_summary.txt" "${cnv_profile_file}"
	mv "${tumor_normal_sample_id}_iCN.seg" "${cnv_segments_file}"
	mv "${tumor_normal_sample_id}_subclonal_cn.txt" "${sclust_subclones_file}"

	Sclust cluster \
	-i "${tumor_normal_sample_id}" \
	-lambda 1e-6 \
	-indel

	mv "${tumor_normal_sample_id}_mclusters.txt" "${mutation_clusters_file}"
	mv "${tumor_normal_sample_id}_mcluster.pdf" "${mutation_clusters_pdf}"
	mv "${tumor_normal_sample_id}_cluster_assignments.txt" "${cluster_assignment_file}"
	"""
}

// END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \\


// ~~~~~~~~~~~~~~~~~ Manta ~~~~~~~~~~~~~~~~~ \\
// START

// Combine all needed reference FASTA files and WGS BED into one channel for use in Manta process
reference_genome_fasta_forManta.combine( reference_genome_fasta_index_forManta )
	.combine( reference_genome_fasta_dict_forManta )
	.combine( gatk_bundle_wgs_bed_forManta )
	.set{ reference_genome_bundle_and_bed_forManta }

// Manta ~ call structural variants and indels from mapped paired-end sequencing reads of matched tumor/normal sample pairs
process svAndIndelCalling_manta {
	publishDir "${params.output_dir}/somatic/manta", mode: 'copy', pattern: '*.{vcf.gz,tbi}'
	tag "${tumor_normal_sample_id}"

	input:
	tuple path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(reference_genome_fasta_forManta), path(reference_genome_fasta_index_forManta), path(reference_genome_fasta_dict_forManta), path(gatk_bundle_wgs_bed_forManta) from tumor_normal_pair_forManta.combine(reference_genome_bundle_and_bed_forManta)

	output:
	tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(candidate_indel_vcf), path(candidate_indel_vcf_index) into bams_and_candidate_indel_vcf_forStrelka
	tuple val(tumor_normal_sample_id), val(tumor_id), path(final_manta_somatic_sv_vcf), path(final_manta_somatic_sv_vcf_index) into manta_sv_vcf_forSurvivorPrep
	tuple path(unfiltered_sv_vcf), path(unfiltered_sv_vcf_index)
	tuple val(tumor_normal_sample_id), path(germline_sv_vcf), path(germline_sv_vcf_index) into germline_indel_vcf_forCaveman

	when:
	params.manta == "on"

	script:
	tumor_id = "${tumor_bam.baseName}".replaceFirst(/\..*$/, "")
	normal_id = "${normal_bam.baseName}".replaceFirst(/\..*$/, "")
	tumor_normal_sample_id = "${tumor_id}_vs_${normal_id}"
	zipped_wgs_bed_forManta = "${gatk_bundle_wgs_bed_forManta}.gz"
	unfiltered_sv_vcf = "${tumor_normal_sample_id}.manta.somatic.sv.unfiltered.vcf.gz"
	unfiltered_sv_vcf_index = "${unfiltered_sv_vcf}.tbi"
	germline_sv_vcf = "${tumor_normal_sample_id}.manta.germline.sv.vcf.gz"
	germline_sv_vcf_index = "${germline_sv_vcf}.tbi"
	candidate_indel_vcf = "${tumor_normal_sample_id}.manta.somatic.indels.unfiltered.vcf.gz"
	candidate_indel_vcf_index = "${candidate_indel_vcf}.tbi"
	final_manta_somatic_sv_vcf = "${tumor_normal_sample_id}.manta.somatic.sv.vcf.gz"
	final_manta_somatic_sv_vcf_index = "${final_manta_somatic_sv_vcf}.tbi"
	"""
	bgzip < "${gatk_bundle_wgs_bed_forManta}" > "${zipped_wgs_bed_forManta}"
	tabix "${zipped_wgs_bed_forManta}"

	python \${MANTA_DIR}/bin/configManta.py \
	--tumorBam "${tumor_bam}" \
	--normalBam "${normal_bam}" \
	--referenceFasta "${reference_genome_fasta_forManta}" \
	--callRegions "${zipped_wgs_bed_forManta}" \
	--runDir manta

	python manta/runWorkflow.py \
	--mode local \
	--jobs ${task.cpus} \
	--memGb ${task.memory.toGiga()}

	mv manta/results/variants/candidateSV.vcf.gz "${unfiltered_sv_vcf}"
	mv manta/results/variants/candidateSV.vcf.gz.tbi "${unfiltered_sv_vcf_index}"
	mv manta/results/variants/candidateSmallIndels.vcf.gz "${candidate_indel_vcf}"
	mv manta/results/variants/candidateSmallIndels.vcf.gz.tbi "${candidate_indel_vcf_index}"

	zcat manta/results/variants/diploidSV.vcf.gz \
	| \
	grep -E "^#|PASS" \
	| \
	bgzip > "${germline_sv_vcf}"
	tabix "${germline_sv_vcf}"

	zcat manta/results/variants/somaticSV.vcf.gz \
	| \
	grep -E "^#|PASS" \
	| \
	bgzip > "${final_manta_somatic_sv_vcf}"
	tabix "${final_manta_somatic_sv_vcf}"
	"""
}

// END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \\


// ~~~~~~~~~~~~~~~~ Strelka2 ~~~~~~~~~~~~~~~ \\
// START

// Combine all needed reference FASTA files and WGS BED into one channel for use in Strelka process
reference_genome_fasta_forStrelka.combine( reference_genome_fasta_index_forStrelka )
	.combine( reference_genome_fasta_dict_forStrelka )
	.combine( gatk_bundle_wgs_bed_forStrelka )
	.set{ reference_genome_bundle_and_bed_forStrelka }

// Strelka2 ~ call germline and somatic small variants from mapped sequencing reads
process snvAndIndelCalling_strelka {
	publishDir "${params.output_dir}/somatic/strelka/rawVcfs", mode: 'symlink', pattern: '*.{vcf.gz,tbi}'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(candidate_indel_vcf), path(candidate_indel_vcf_index), path(reference_genome_fasta_forStrelka), path(reference_genome_fasta_index_forStrelka), path(reference_genome_fasta_dict_forStrelka), path(gatk_bundle_wgs_bed_forStrelka) from bams_and_candidate_indel_vcf_forStrelka.combine(reference_genome_bundle_and_bed_forStrelka)

	output:
	tuple val(tumor_normal_sample_id), path(unfiltered_strelka_snv_vcf), path(unfiltered_strelka_snv_vcf_index), path(unfiltered_strelka_indel_vcf), path(unfiltered_strelka_indel_vcf_index) into unfiltered_snv_and_indel_vcfs_forStrelkaBcftools

	when:
	params.strelka == "on" && params.manta == "on"

	script:
	zipped_wgs_bed_forStrelka = "${gatk_bundle_wgs_bed_forStrelka}.gz"
	unfiltered_strelka_snv_vcf = "${tumor_normal_sample_id}.strelka.somatic.snv.unfiltered.vcf.gz"
	unfiltered_strelka_snv_vcf_index = "${unfiltered_strelka_snv_vcf}.tbi"
	unfiltered_strelka_indel_vcf = "${tumor_normal_sample_id}.strelka.somatic.indel.unfiltered.vcf.gz"
	unfiltered_strelka_indel_vcf_index = "${unfiltered_strelka_indel_vcf}.tbi"
	"""
	bgzip < "${gatk_bundle_wgs_bed_forStrelka}" > "${zipped_wgs_bed_forStrelka}"
	tabix "${zipped_wgs_bed_forStrelka}"

	python \${STRELKA_DIR}/bin/configureStrelkaSomaticWorkflow.py \
	--normalBam "${normal_bam}" \
	--tumorBam "${tumor_bam}" \
	--referenceFasta "${reference_genome_fasta_forStrelka}" \
	--indelCandidates "${candidate_indel_vcf}" \
	--callRegions "${zipped_wgs_bed_forStrelka}" \
	--runDir strelka

	python strelka/runWorkflow.py \
	--mode local \
	--jobs ${task.cpus} \
	--memGb ${task.memory.toGiga()}

	mv strelka/results/variants/somatic.snvs.vcf.gz "${unfiltered_strelka_snv_vcf}"
	mv strelka/results/variants/somatic.snvs.vcf.gz.tbi "${unfiltered_strelka_snv_vcf_index}"
	mv strelka/results/variants/somatic.indels.vcf.gz "${unfiltered_strelka_indel_vcf}"
	mv strelka/results/variants/somatic.indels.vcf.gz.tbi "${unfiltered_strelka_indel_vcf_index}"
	"""
}

// Combine all needed reference FASTA files into one channel for use in Strelka / BCFtools Norm process
reference_genome_fasta_forStrelkaBcftools.combine( reference_genome_fasta_index_forStrelkaBcftools )
	.combine( reference_genome_fasta_dict_forStrelkaBcftools )
	.set{ reference_genome_bundle_forStrelkaBcftools }

// BCFtools norm ~ split multiallelic sites into multiple rows then left-align and normalize indels
process splitMultiallelicAndLeftNormalizeStrelkaVcf_bcftools {
	publishDir "${params.output_dir}/somatic/strelka", mode: 'copy', pattern: '*.{vcf.gz,tbi,txt}'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(unfiltered_strelka_snv_vcf), path(unfiltered_strelka_snv_vcf_index), path(unfiltered_strelka_indel_vcf), path(unfiltered_strelka_indel_vcf_index), path(reference_genome_fasta_forStrelkaBcftools), path(reference_genome_fasta_index_forStrelkaBcftools), path(reference_genome_fasta_dict_forStrelkaBcftools) from unfiltered_snv_and_indel_vcfs_forStrelkaBcftools.combine(reference_genome_bundle_forStrelkaBcftools)

	output:
	tuple val(tumor_normal_sample_id), path(final_strelka_snv_vcf), path(final_strelka_snv_vcf_index) into final_strelka_snv_vcf_forConsensus
	tuple val(tumor_normal_sample_id), path(final_strelka_indel_vcf), path(final_strelka_indel_vcf_index) into final_strelka_indel_vcf_forConsensus
	path strelka_snv_multiallelics_stats
	path strelka_indel_multiallelics_stats
	path strelka_indel_realign_normalize_stats

	when:
	params.strelka == "on" && params.manta == "on"

	script:
	final_strelka_snv_vcf = "${tumor_normal_sample_id}.strelka.somatic.snv.vcf.gz"
	final_strelka_snv_vcf_index ="${final_strelka_snv_vcf}.tbi"
	final_strelka_indel_vcf = "${tumor_normal_sample_id}.strelka.somatic.indel.vcf.gz"
	final_strelka_indel_vcf_index = "${final_strelka_indel_vcf}.tbi"
	strelka_snv_multiallelics_stats = "${tumor_normal_sample_id}.strelka.snv.multiallelicsstats.txt"
	strelka_indel_multiallelics_stats = "${tumor_normal_sample_id}.strelka.indel.multiallelicsstats.txt"
	strelka_indel_realign_normalize_stats = "${tumor_normal_sample_id}.strelka.indel.realignnormalizestats.txt"
	"""
	zgrep -E "^#|PASS" "${unfiltered_strelka_snv_vcf}" \
	| \
	bcftools norm \
	--threads ${task.cpus} \
	--multiallelics -snps \
	--output-type z \
	--output "${final_strelka_snv_vcf}" \
	- 2>"${strelka_snv_multiallelics_stats}"

	tabix "${final_strelka_snv_vcf}"
	
	zgrep -E "^#|PASS" "${unfiltered_strelka_indel_vcf}" \
	| \
	bcftools norm \
	--threads ${task.cpus} \
	--multiallelics -indels \
	--output-type z \
	- 2>"${strelka_indel_multiallelics_stats}" \
	| \
	bcftools norm \
	--threads ${task.cpus} \
	--fasta-ref "${reference_genome_fasta_forStrelkaBcftools}" \
	--output-type z \
	--output "${final_strelka_indel_vcf}" \
	- 2>"${strelka_indel_realign_normalize_stats}"

	tabix "${final_strelka_indel_vcf}"
	"""
}

// END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \\


// ~~~~~~~~~~~~~~~~ CaVEMan ~~~~~~~~~~~~~~~ \\
// START

// BEDOPS vcf2bed / sort-bed ~ convert germline indel VCF to sorted BED file
process prepGermlineBedForCaveman_bedops {
	publishDir "${params.output_dir}/somatic/caveman/intermediates", mode: 'symlink', pattern: '*.{bed,tbi}'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(germline_sv_vcf), path(germline_sv_vcf_index) from germline_indel_vcf_forCaveman

	output:
	tuple val(tumor_normal_sample_id), path(germline_indel_bed), path(germline_indel_bed_index) into germline_indel_bed_forCaveman

	when:
	params.caveman == "on" && params.ascatngs == "on" && params.manta == "on"

	script:
	germline_insertions_bed = "${tumor_normal_sample_id}.caveman.germline.ins.bed"
	germline_deletions_bed = "${tumor_normal_sample_id}.caveman.germline.dels.bed"
	germline_indel_bed = "${tumor_normal_sample_id}.caveman.germline.indels.bed.gz"
	germline_indel_bed_index = "${germline_indel_bed}.tbi"
	"""
	zgrep -E '^#|MantaINS.*PASS' "${germline_sv_vcf}" \
	| \
	vcf2bed --insertions \
	| \
	cut -f 1-4 > "${germline_insertions_bed}"

	zgrep -E '^#|MantaDEL.*PASS' "${germline_sv_vcf}" \
	| \
	vcf2bed --deletions \
	| \
	cut -f 1-4 > "${germline_deletions_bed}"

	sort-bed "${germline_insertions_bed}" "${germline_deletions_bed}" \
	| \
	bgzip > "${germline_indel_bed}"

	tabix -p bed "${germline_indel_bed}"
	"""
}

// Combine all needed reference FASTA files, blacklist BED, and flagging resource BEDs into one channel for CaVEMan processes
reference_genome_fasta_forCaveman.combine( reference_genome_fasta_index_forCaveman )
	.combine( reference_genome_fasta_dict_forCaveman )
	.combine( gatk_bundle_wgs_bed_blacklist_1based_forCaveman )
	.combine( centromeric_repeats_bed )
	.combine( centromeric_repeats_bed_index )
	.combine( simple_repeats_bed )
	.combine( simple_repeats_bed_index )
	.combine( dbsnp_bed )
	.combine( dbsnp_bed_index )
	.combine( unmatched_normal_bed )
	.combine( unmatched_normal_bed_index )
	.set{ resource_bundle_forCaveman }

// CaVEMan setup ~ generate configuration files for all execution steps
process setup_caveman {
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(cnv_profile_final), path(run_statistics), path(germline_indel_bed), path(germline_indel_bed_index), path(reference_genome_fasta_forCaveman), path(reference_genome_fasta_index_forCaveman), path(reference_genome_fasta_dict_forCaveman), path(gatk_bundle_wgs_bed_blacklist_1based_forCaveman), path(unmatched_normal_bed), path(unmatched_normal_bed_index), path(centromeric_repeats_bed), path(centromeric_repeats_bed_index), path(simple_repeats_bed), path(simple_repeats_bed_index), path(dbsnp_bed), path(dbsnp_bed_index) from bams_cnv_profile_and_statistics_forCaveman.join(germline_indel_bed_forCaveman).combine(resource_bundle_forCaveman)

	output:
	tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(tumor_cnv_profile_bed), path(normal_cnv_profile_bed), path(run_statistics), path(germline_indel_bed), path(germline_indel_bed_index), path(reference_genome_fasta_forCaveman), path(reference_genome_fasta_index_forCaveman), path(reference_genome_fasta_dict_forCaveman), path(gatk_bundle_wgs_bed_blacklist_1based_forCaveman), path(unmatched_normal_bed), path(unmatched_normal_bed_index), path(centromeric_repeats_bed), path(centromeric_repeats_bed_index), path(simple_repeats_bed), path(simple_repeats_bed_index), path(dbsnp_bed), path(dbsnp_bed_index), path(postprocessing_config_file), path(config_file), path(alg_bean_file) into setup_forCavemanSplit, setup_forCavemanConcat, setup_forCavemanMstep

	when:
	params.caveman == "on" && params.ascatngs == "on" && params.manta == "on"

	script:
	tumor_cnv_profile_bed = "${tumor_normal_sample_id}.caveman.tumor.cvn.bed"
	normal_cnv_profile_bed = "${tumor_normal_sample_id}.caveman.normal.cvn.bed"
	postprocessing_config_file = "${tumor_normal_sample_id}.cavemanpostprocessing.config.ini"
	config_species = "HOMO_SAPIENS"
	config_study_type = "WGS"
	config_file = "tmpCaveman/caveman.cfg.ini"
	alg_bean_file = "tmpCaveman/alg_bean"
	"""
	if [[ ! "${tumor_bam_index}" =~ .bam.bai\$ ]]; then
		cp "${tumor_bam_index}" "${tumor_bam}.bai"
	fi
	if [[ ! "${normal_bam_index}" =~ .bam.bai\$ ]]; then
		cp "${normal_bam_index}" "${normal_bam}.bai"
	fi

	touch "${postprocessing_config_file}"
	echo "[${config_species}_${config_study_type} PARAMS]" >> "${postprocessing_config_file}"
	echo "keepSW=1" >> "${postprocessing_config_file}"
	echo "minAnalysedQual=11" >> "${postprocessing_config_file}"
	echo "maxMatchedNormalAlleleProportion=0.05" >> "${postprocessing_config_file}"
	echo "maxPhasingMinorityStrandReadProportion=0.04" >> "${postprocessing_config_file}"
	echo "readPosBeginningOfReadIgnoreProportion=0.08" >> "${postprocessing_config_file}"
	echo "readPosTwoThirdsOfReadExtendProportion=0.08" >> "${postprocessing_config_file}"
	echo "matchedNormalAlleleHiCvgCutoff=2" >> "${postprocessing_config_file}"
	echo "maxMatchedNormalAlleleHiCvgProportion=0.03" >> "${postprocessing_config_file}"
	echo "pentamerMinPassAvgQual=20" >> "${postprocessing_config_file}"
	echo "samePosMaxPercent=80" >> "${postprocessing_config_file}"
	echo "maxTumIndelProportion=10" >> "${postprocessing_config_file}"
	echo "maxNormIndelProportion=10" >> "${postprocessing_config_file}"
	echo "minPassAvgMapQual=21" >> "${postprocessing_config_file}"
	echo "minPassPhaseQual=21" >> "${postprocessing_config_file}"
	echo "minDepthQual=25" >> "${postprocessing_config_file}"
	echo "minNormMutAllelequal=15" >> "${postprocessing_config_file}"
	echo "minRdPosDepth=8" >> "${postprocessing_config_file}"
	echo "vcfUnmatchedMinMutAlleleCvg=3" >> "${postprocessing_config_file}"
	echo "vcfUnmatchedMinSamplePct=1" >> "${postprocessing_config_file}"
	echo "matchedNormalMaxMutProportion=0.2" >> "${postprocessing_config_file}"
	echo "minSingleEndCoverage=10" >> "${postprocessing_config_file}"
	echo "" >> "${postprocessing_config_file}"

	echo "[${config_species}_${config_study_type} FLAGLIST]" >> "${postprocessing_config_file}"
	echo "flagList=<<LST" >> "${postprocessing_config_file}"
	echo "depthFlag" >> "${postprocessing_config_file}"
	echo "readPositionFlag" >> "${postprocessing_config_file}"
	echo "matchedNormalFlag" >> "${postprocessing_config_file}"
	echo "pentamericMotifFlag" >> "${postprocessing_config_file}"
	echo "avgMapQualFlag" >> "${postprocessing_config_file}"
	echo "simpleRepeatFlag" >> "${postprocessing_config_file}"
	echo "centromericRepeatFlag" >> "${postprocessing_config_file}"
	echo "snpFlag" >> "${postprocessing_config_file}"
	echo "phasingFlag" >> "${postprocessing_config_file}"
	echo "germlineIndelFlag" >> "${postprocessing_config_file}"
	echo "unmatchedNormalVcfFlag" >> "${postprocessing_config_file}"
	echo "singleEndFlag" >> "${postprocessing_config_file}"
	echo "matchedNormalProportion" >> "${postprocessing_config_file}"
	echo "alignmentScoreReadLengthAdjustedFlag" >> "${postprocessing_config_file}"
	echo "clippingMedianFlag" >> "${postprocessing_config_file}"
	echo "alnScoreMedianFlag" >> "${postprocessing_config_file}"
	echo "LST" >> "${postprocessing_config_file}"
	echo "" >> "${postprocessing_config_file}"

	echo "[${config_species}_${config_study_type} BEDFILES]" >> "${postprocessing_config_file}"
	echo "centromericRepeatBed=${centromeric_repeats_bed}" >> "${postprocessing_config_file}"
	echo "simpleRepeatBed=${simple_repeats_bed}" >> "${postprocessing_config_file}"
	echo "snpBed=${dbsnp_bed}" >> "${postprocessing_config_file}"
	echo "" >> "${postprocessing_config_file}"

	grep -v 'segment_number' "${cnv_profile_final}" \
	| \
	awk -F, 'BEGIN {OFS="\t"} {print \$2,\$3,\$4,\$7}' > "${tumor_cnv_profile_bed}"

	grep -v 'segment_number' "${cnv_profile_final}" \
	| \
	awk -F, 'BEGIN {OFS="\t"} {print \$2,\$3,\$4,\$5}' > "${normal_cnv_profile_bed}"

	caveman.pl \
	-outdir . \
	-reference "${reference_genome_fasta_index_forCaveman}" \
	-tumour-bam "${tumor_bam}" \
	-normal-bam "${normal_bam}" \
	-ignore-file "${gatk_bundle_wgs_bed_blacklist_1based_forCaveman}" \
	-tumour-cn "${tumor_cnv_profile_bed}" \
	-normal-cn "${normal_cnv_profile_bed}" \
	-species Homo_sapiens \
	-species-assembly GRCh38 \
	-flag-bed-files . \
	-germline-indel "${germline_indel_bed}" \
	-unmatched-vcf "${unmatched_normal_bed}" \
	-seqType genome \
	-threads ${task.cpus} \
	-normal-contamination "${run_statistics}" \
	-flagConfig "${postprocessing_config_file}" \
	-process setup \
	-index 1
	"""
}

// CaVEMan split ~ split the genome into chunks by readsize and hard stop forced by contig ends
process split_caveman {
	tag "C=${chromosome} ${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(tumor_cnv_profile_bed), path(normal_cnv_profile_bed), path(run_statistics), path(germline_indel_bed), path(germline_indel_bed_index), path(reference_genome_fasta_forCaveman), path(reference_genome_fasta_index_forCaveman), path(reference_genome_fasta_dict_forCaveman), path(gatk_bundle_wgs_bed_blacklist_1based_forCaveman), path(unmatched_normal_bed), path(unmatched_normal_bed_index), path(centromeric_repeats_bed), path(centromeric_repeats_bed_index), path(simple_repeats_bed), path(simple_repeats_bed_index), path(dbsnp_bed), path(dbsnp_bed_index), path(postprocessing_config_file), path(config_file), path(alg_bean_file) from setup_forCavemanSplit
	each chromosome from chromosome_list_forCavemanSplit

	output:
	tuple val(tumor_normal_sample_id), path(split_list_per_chromosome), path(read_position_per_chromosome) into split_per_chromosome_forCavemanMstep

	//tuple val(tumor_normal_sample_id), path(split_list_per_chromosome), path(read_position_per_chromosome) into split_per_chromosome_forCavemanConcat, split_per_chromosome_forCavemanMstep

	when:
	params.caveman == "on" && params.ascatngs == "on" && params.manta == "on"

	script:
	split_list_per_chromosome = "tmpCaveman/splitList.${chromosome}"
	read_position_per_chromosome = "tmpCaveman/readpos.${chromosome}"
	"""
	if [[ ! "${tumor_bam_index}" =~ .bam.bai\$ ]]; then
		cp "${tumor_bam_index}" "${tumor_bam}.bai"
	fi
	if [[ ! "${normal_bam_index}" =~ .bam.bai\$ ]]; then
		cp "${normal_bam_index}" "${normal_bam}.bai"
	fi

	sed -i'' 's|CWD=.*|CWD='"\$PWD"'|' "${config_file}"
	sed -i'' 's|ALG_FILE=.*|ALG_FILE='"\$PWD/${alg_bean_file}"'|' "${config_file}"

	mkdir -p tmpCaveman/
	mv "${config_file}" tmpCaveman/

	i=\$(grep -wn "${chromosome}" "${reference_genome_fasta_index_forCaveman}" | cut -f 1 | cut -d ':' -f 1)

	caveman.pl \
	-outdir . \
	-reference "${reference_genome_fasta_index_forCaveman}" \
	-tumour-bam "${tumor_bam}" \
	-normal-bam "${normal_bam}" \
	-ignore-file "${gatk_bundle_wgs_bed_blacklist_1based_forCaveman}" \
	-tumour-cn "${tumor_cnv_profile_bed}" \
	-normal-cn "${normal_cnv_profile_bed}" \
	-species Homo_sapiens \
	-species-assembly GRCh38 \
	-flag-bed-files . \
	-germline-indel "${germline_indel_bed}" \
	-unmatched-vcf "${unmatched_normal_bed}" \
	-seqType genome \
	-threads ${task.cpus} \
	-normal-contamination "${run_statistics}" \
	-flagConfig "${postprocessing_config_file}" \
	-process split \
	-index \${i}
	"""
}

// CaVEMan mstep ~ build a profile of each split section of the genome using various covariates
process mstep_caveman {
	tag "C=${chromosome} ${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(tumor_cnv_profile_bed), path(normal_cnv_profile_bed), path(run_statistics), path(germline_indel_bed), path(germline_indel_bed_index), path(reference_genome_fasta_forCaveman), path(reference_genome_fasta_index_forCaveman), path(reference_genome_fasta_dict_forCaveman), path(gatk_bundle_wgs_bed_blacklist_1based_forCaveman), path(unmatched_normal_bed), path(unmatched_normal_bed_index), path(centromeric_repeats_bed), path(centromeric_repeats_bed_index), path(simple_repeats_bed), path(simple_repeats_bed_index), path(dbsnp_bed), path(dbsnp_bed_index), path(postprocessing_config_file), path(config_file), path(alg_bean_file), path(split_list_per_chromosome), path(read_position_per_chromosome) from setup_forCavemanMstep.join(split_per_chromosome_forCavemanMstep.groupTuple())
	each chromosome from chromosome_list_forCavemanMstep

	//output:


	when:
	params.caveman == "on" && params.ascatngs == "on" && params.manta == "on"

	script:
	"""
	if [[ ! "${tumor_bam_index}" =~ .bam.bai\$ ]]; then
		cp "${tumor_bam_index}" "${tumor_bam}.bai"
	fi
	if [[ ! "${normal_bam_index}" =~ .bam.bai\$ ]]; then
		cp "${normal_bam_index}" "${normal_bam}.bai"
	fi

	sed -i'' 's|CWD=.*|CWD='"\$PWD"'|' "${config_file}"
	sed -i'' 's|ALG_FILE=.*|ALG_FILE='"\$PWD/${alg_bean_file}"'|' "${config_file}"

	mkdir -p tmpCaveman/
	mv "${config_file}" tmpCaveman/
	cp "splitList.${chromosome}" tmpCaveman/splitList
	mv "readpos.${chromosome}" tmpCaveman/

	caveman.pl \
	-outdir . \
	-reference "${reference_genome_fasta_index_forCaveman}" \
	-tumour-bam "${tumor_bam}" \
	-normal-bam "${normal_bam}" \
	-ignore-file "${gatk_bundle_wgs_bed_blacklist_1based_forCaveman}" \
	-tumour-cn "${tumor_cnv_profile_bed}" \
	-normal-cn "${normal_cnv_profile_bed}" \
	-species Homo_sapiens \
	-species-assembly GRCh38 \
	-flag-bed-files . \
	-germline-indel "${germline_indel_bed}" \
	-unmatched-vcf "${unmatched_normal_bed}" \
	-seqType genome \
	-threads ${task.cpus} \
	-normal-contamination "${run_statistics}" \
	-flagConfig "${postprocessing_config_file}" \
	-process mstep
	done
	"""
}

// END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \\


// ~~~~~~~~~~~~~~~~~ SvABA ~~~~~~~~~~~~~~~~~ \\
// START

// Combine all needed reference FASTA files, WGS BED, and reference VCF files into one channel for use in SvABA process
bwa_ref_genome_files.collect()
	.combine( reference_genome_fasta_dict_forSvaba )
	.combine( gatk_bundle_wgs_bed_forSvaba )
	.set{ bwa_ref_genome_and_wgs_bed }

bwa_ref_genome_and_wgs_bed.combine( dbsnp_known_indel_ref_vcf )
	.combine( dbsnp_known_indel_ref_vcf_index )
	.set{ bwa_ref_genome_wgs_bed_and_ref_vcf }

// SvABA ~ detecting structural variants using genome-wide local assembly
process svAndIndelCalling_svaba {
	publishDir "${params.output_dir}/somatic/svaba", mode: 'copy', pattern: '*.{somatic.sv.vcf.gz,somatic.sv.vcf.gz.tbi}'
	publishDir "${params.output_dir}/somatic/svaba/intermediates", mode: 'symlink', pattern: '*.{germline,unfiltered}*'
	tag "${tumor_normal_sample_id}"

	input:
	tuple path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path("Homo_sapiens_assembly38.fasta"), path("Homo_sapiens_assembly38.fasta.fai"), path("Homo_sapiens_assembly38.fasta.64.alt"), path("Homo_sapiens_assembly38.fasta.64.amb"), path("Homo_sapiens_assembly38.fasta.64.ann"), path("Homo_sapiens_assembly38.fasta.64.bwt"), path("Homo_sapiens_assembly38.fasta.64.pac"), path("Homo_sapiens_assembly38.fasta.64.sa"), path(reference_genome_fasta_dict_forSvaba), path(gatk_bundle_wgs_bed_forSvaba), path(dbsnp_known_indel_ref_vcf), path(dbsnp_known_indel_ref_vcf_index) from tumor_normal_pair_forSvaba.combine(bwa_ref_genome_wgs_bed_and_ref_vcf)

	output:
	tuple val(tumor_normal_sample_id), path(filtered_somatic_indel_vcf), path(filtered_somatic_indel_vcf_index) into filtered_indel_vcf_forSvabaBcftools
	tuple val(tumor_normal_sample_id), val(tumor_id), path(final_svaba_somatic_sv_vcf), path(final_svaba_somatic_sv_vcf_index), path(sample_renaming_file) into svaba_sv_vcf_forSurvivorPrep
	path unfiltered_somatic_indel_vcf
	path unfiltered_somatic_sv_vcf
	path germline_indel_vcf
	path germline_sv_vcf

	when:
	params.svaba == "on"

	script:
	tumor_id = "${tumor_bam.baseName}".replaceFirst(/\..*$/, "")
	normal_id = "${normal_bam.baseName}".replaceFirst(/\..*$/, "")
	tumor_normal_sample_id = "${tumor_id}_vs_${normal_id}"
	contig_alignment_plot = "${tumor_normal_sample_id}.svaba.alignments.txt.gz"
	germline_indel_vcf = "${tumor_normal_sample_id}.svaba.germline.indel.vcf.gz"
	germline_sv_vcf = "${tumor_normal_sample_id}.svaba.germline.sv.vcf.gz"
	unfiltered_somatic_indel_vcf = "${tumor_normal_sample_id}.svaba.somatic.unfiltered.indel.vcf.gz"
	unfiltered_somatic_sv_vcf = "${tumor_normal_sample_id}.svaba.somatic.unfiltered.sv.vcf.gz"
	filtered_somatic_indel_vcf = "${tumor_normal_sample_id}.svaba.somatic.filtered.indel.vcf.gz"
	filtered_somatic_indel_vcf_index = "${filtered_somatic_indel_vcf}.tbi"
	final_svaba_somatic_sv_vcf = "${tumor_normal_sample_id}.svaba.somatic.sv.vcf.gz"
	final_svaba_somatic_sv_vcf_index = "${final_svaba_somatic_sv_vcf}.tbi"
	sample_renaming_file = "sample_renaming_file.txt"
	"""
	svaba run \
	-t "${tumor_bam}" \
	-n "${normal_bam}" \
	--reference-genome Homo_sapiens_assembly38.fasta \
	--region "${gatk_bundle_wgs_bed_forSvaba}" \
	--id-string "${tumor_normal_sample_id}" \
	--dbsnp-vcf "${dbsnp_known_indel_ref_vcf}" \
	--threads "${task.cpus}" \
	--verbose 1 \
	--g-zip

	mv "${tumor_normal_sample_id}.alignments.txt.gz" "${contig_alignment_plot}"
	mv "${tumor_normal_sample_id}.svaba.unfiltered.somatic.indel.vcf.gz" "${unfiltered_somatic_indel_vcf}"
	mv "${tumor_normal_sample_id}.svaba.unfiltered.somatic.sv.vcf.gz" "${unfiltered_somatic_sv_vcf}"
	mv "${tumor_normal_sample_id}.svaba.somatic.indel.vcf.gz" "${filtered_somatic_indel_vcf}"

	tabix "${filtered_somatic_indel_vcf}"
	tabix "${final_svaba_somatic_sv_vcf}"

	touch "${sample_renaming_file}"
	echo "${tumor_bam} ${tumor_id}" >> "${sample_renaming_file}"
	"""
}

// Combine all needed reference FASTA files into one channel for use in SvABA / BCFtools Norm process
reference_genome_fasta_forSvabaBcftools.combine( reference_genome_fasta_index_forSvabaBcftools )
	.combine( reference_genome_fasta_dict_forSvabaBcftools )
	.set{ reference_genome_bundle_forSvabaBcftools }

// BCFtools Norm ~ left-align and normalize indels
process leftNormalizeSvabaVcf_bcftools {
	publishDir "${params.output_dir}/somatic/svaba", mode: 'copy', pattern: '*.{vcf.gz,tbi,txt}'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(filtered_somatic_indel_vcf), path(filtered_somatic_indel_vcf_index), path(reference_genome_fasta_forSvabaBcftools), path(reference_genome_fasta_index_forSvabaBcftools),
	path(reference_genome_fasta_dict_forSvabaBcftools) from filtered_indel_vcf_forSvabaBcftools.combine(reference_genome_bundle_forSvabaBcftools)

	output:
	tuple val(tumor_normal_sample_id), path(final_svaba_indel_vcf), path(final_svaba_indel_vcf_index) into final_svaba_indel_vcf_forConsensus
	path svaba_realign_normalize_stats

	when:
	params.svaba == "on"

	script:
	final_svaba_indel_vcf = "${tumor_normal_sample_id}.svaba.somatic.indel.vcf.gz"
	final_svaba_indel_vcf_index = "${final_svaba_indel_vcf}.tbi"
	svaba_realign_normalize_stats = "${tumor_normal_sample_id}.svaba.realignnormalizestats.txt"
	"""
	bcftools norm \
	--threads ${task.cpus} \
	--fasta-ref "${reference_genome_fasta_forSvabaBcftools}" \
	--output-type z \
	--output "${final_svaba_indel_vcf}" \
	"${filtered_somatic_indel_vcf}" 2>"${svaba_realign_normalize_stats}"

	tabix "${final_svaba_indel_vcf}"
	"""
}

// END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \\


// ~~~~~~~~~~~~~~~~ DELLY2 ~~~~~~~~~~~~~~~~~ \\
// START

// Combine all reference FASTA files and calling blacklist into one channel for use in DELLY2 process
reference_genome_fasta_forDelly.combine( reference_genome_fasta_index_forDelly )
	.combine( reference_genome_fasta_dict_forDelly )
	.combine( gatk_bundle_wgs_bed_blacklist_0based_forDelly )
	.set{ reference_genome_and_blacklist_bundle_forDelly }

// DELLY2 ~ discover structural variants using paired-ends, split-reads and read-depth
process svAndIndelCalling_delly {
	publishDir "${params.output_dir}/somatic/delly", mode: 'copy', pattern: '*.{vcf.gz,tbi}'
	tag "${tumor_normal_sample_id}"

	input:
	tuple path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(reference_genome_fasta_forDelly), path(reference_genome_fasta_index_forDelly), path(reference_genome_fasta_dict_forDelly), path(gatk_bundle_wgs_bed_blacklist_0based_forDelly) from tumor_normal_pair_forDelly.combine(reference_genome_and_blacklist_bundle_forDelly)

	output:
	tuple val(tumor_normal_sample_id), val(tumor_id), path(final_delly_somatic_sv_vcf), path(final_delly_somatic_sv_vcf_index) into delly_sv_vcf_forSurvivorPrep

	when:
	params.delly == "on"

	script:
	tumor_id = "${tumor_bam.baseName}".replaceFirst(/\..*$/, "")
	normal_id = "${normal_bam.baseName}".replaceFirst(/\..*$/, "")
	tumor_normal_sample_id = "${tumor_id}_vs_${normal_id}"
	final_delly_somatic_sv_vcf = "${tumor_normal_sample_id}.delly.somatic.sv.vcf.gz"
	final_delly_somatic_sv_vcf_index = "${final_delly_somatic_sv_vcf}.tbi"
	"""
	delly call \
	--genome "${reference_genome_fasta_forDelly}" \
	--exclude "${gatk_bundle_wgs_bed_blacklist_0based_forDelly}" \
	--outfile "${tumor_normal_sample_id}.delly.somatic.sv.unfiltered.bcf" \
	"${tumor_bam}" "${normal_bam}"

	touch samples.tsv
	echo "${tumor_id}\ttumor" >> samples.tsv
	echo "${normal_id}\tcontrol" >> samples.tsv

	delly filter \
	--filter somatic \
	--pass \
	--samples samples.tsv \
	--outfile "${tumor_normal_sample_id}.delly.somatic.sv.bcf" \
	"${tumor_normal_sample_id}.delly.somatic.sv.unfiltered.bcf"

	bcftools view \
	--output-type z \
	"${tumor_normal_sample_id}.delly.somatic.sv.bcf" > "${final_delly_somatic_sv_vcf}"

	tabix "${final_delly_somatic_sv_vcf}"
	"""
}

// END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \\


// ~~~~~~~~~~~~~~~~ GRIDSS ~~~~~~~~~~~~~~~~~ \\
// START

// Combine all needed reference FASTA files and WGS blacklist BED into one channel for GRIDSS setupreference process
reference_genome_fasta_forGridssSetup.combine( reference_genome_fasta_index_forGridssSetup )
	.combine( gatk_bundle_wgs_bed_blacklist_1based_forGridssSetup )
	.set{ reference_genome_bundle_and_bed_forGridssSetup }

// GRIDSS setupreference ~ generate additional needed files in the same directory as the reference genome
process setupreference_gridss {

	input:
	tuple path(reference_genome_fasta_forGridssSetup), path(reference_genome_fasta_index_forGridssSetup), path(gatk_bundle_wgs_bed_blacklist_1based_forGridssSetup) from reference_genome_bundle_and_bed_forGridssSetup

	output:
	tuple path(reference_genome_fasta_forGridssSetup), path(reference_genome_fasta_index_forGridssSetup), path(gatk_bundle_wgs_bed_blacklist_1based_forGridssSetup), path("Homo_sapiens_assembly38.fasta.amb"), path("Homo_sapiens_assembly38.fasta.ann"), path("Homo_sapiens_assembly38.fasta.bwt"), path("Homo_sapiens_assembly38.fasta.dict"), path("Homo_sapiens_assembly38.fasta.gridsscache"), path("Homo_sapiens_assembly38.fasta.img"), path("Homo_sapiens_assembly38.fasta.pac"), path("Homo_sapiens_assembly38.fasta.sa") into setup_reference_output_forGridssPreprocess

	when:
	params.gridss == "on"

	script:
	"""
	gridss.sh --steps setupreference \
	--jvmheap "${task.memory.toGiga()}"g \
	--otherjvmheap "${task.memory.toGiga()}"g \
	--threads "${task.cpus}" \
	--reference "${reference_genome_fasta_forGridssSetup}" \
	--blacklist "${gatk_bundle_wgs_bed_blacklist_1based_forGridssSetup}"
	"""
}

// GRIDSS preprocess ~ preprocess input tumor and normal BAMs separately
process preprocess_gridss {
	publishDir "${params.output_dir}/somatic/gridss/intermediates", mode: 'symlink', pattern: '*${tumor_bam_working_dir}'
	publishDir "${params.output_dir}/somatic/gridss/intermediates", mode: 'symlink', pattern: '*${normal_bam_working_dir}'
	tag "T=${tumor_id} N=${normal_id}"

	input:
	tuple path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(reference_genome_fasta_forGridssSetup), path(reference_genome_fasta_index_forGridssSetup), path(gatk_bundle_wgs_bed_blacklist_1based_forGridssSetup), path("Homo_sapiens_assembly38.fasta.amb"), path("Homo_sapiens_assembly38.fasta.ann"), path("Homo_sapiens_assembly38.fasta.bwt"), path("Homo_sapiens_assembly38.fasta.dict"), path("Homo_sapiens_assembly38.fasta.gridsscache"), path("Homo_sapiens_assembly38.fasta.img"), path("Homo_sapiens_assembly38.fasta.pac"), path("Homo_sapiens_assembly38.fasta.sa") from tumor_normal_pair_forGridssPreprocess.combine(setup_reference_output_forGridssPreprocess)

	output:
	tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(reference_genome_fasta_forGridssSetup), path(reference_genome_fasta_index_forGridssSetup), path(gatk_bundle_wgs_bed_blacklist_1based_forGridssSetup), path("Homo_sapiens_assembly38.fasta.amb"), path("Homo_sapiens_assembly38.fasta.ann"), path("Homo_sapiens_assembly38.fasta.bwt"), path("Homo_sapiens_assembly38.fasta.dict"), path("Homo_sapiens_assembly38.fasta.gridsscache"), path("Homo_sapiens_assembly38.fasta.img"), path("Homo_sapiens_assembly38.fasta.pac"), path("Homo_sapiens_assembly38.fasta.sa"), path(tumor_bam_working_dir), path(normal_bam_working_dir) into preprocess_output_forGridssAssemble, preprocess_output_forGridssMergeAssembly

	when:
	params.gridss == "on"

	script:
	tumor_id = "${tumor_bam.baseName}".replaceFirst(/\..*$/, "")
	normal_id = "${normal_bam.baseName}".replaceFirst(/\..*$/, "")
	tumor_normal_sample_id = "${tumor_id}_vs_${normal_id}"
	tumor_bam_working_dir = "${tumor_bam}.gridss.working"
	normal_bam_working_dir = "${normal_bam}.gridss.working"
	"""
	gridss.sh --steps preprocess \
	--jvmheap "${task.memory.toGiga()}"g \
	--otherjvmheap "${task.memory.toGiga()}"g \
	--threads "${task.cpus}" \
	--reference "${reference_genome_fasta_forGridssSetup}" \
	--blacklist "${gatk_bundle_wgs_bed_blacklist_1based_forGridssSetup}" \
	"${tumor_bam}"

	gridss.sh --steps preprocess \
	--jvmheap "${task.memory.toGiga()}"g \
	--otherjvmheap "${task.memory.toGiga()}"g \
	--threads "${task.cpus}" \
	--reference "${reference_genome_fasta_forGridssSetup}" \
	--blacklist "${gatk_bundle_wgs_bed_blacklist_1based_forGridssSetup}" \
	"${normal_bam}"
	"""
}

// Create channel for job index of each GRIDSS assemble processs
Channel
	.from( 0,1 )
	.set{ job_index_list_forGridssAssemble }

// GRIDSS assemble ~ perform breakend assembly split up across multiple jobs
process assemble_gridss {
	publishDir "${params.output_dir}/somatic/gridss/intermediates", mode: 'symlink', pattern: '*${assembly_bam_per_index_working_dir}'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(reference_genome_fasta_forGridssSetup), path(reference_genome_fasta_index_forGridssSetup), path(gatk_bundle_wgs_bed_blacklist_1based_forGridssSetup), path("Homo_sapiens_assembly38.fasta.amb"), path("Homo_sapiens_assembly38.fasta.ann"), path("Homo_sapiens_assembly38.fasta.bwt"), path("Homo_sapiens_assembly38.fasta.dict"), path("Homo_sapiens_assembly38.fasta.gridsscache"), path("Homo_sapiens_assembly38.fasta.img"), path("Homo_sapiens_assembly38.fasta.pac"), path("Homo_sapiens_assembly38.fasta.sa"), path(tumor_bam_working_dir), path(normal_bam_working_dir) from preprocess_output_forGridssAssemble
	each index from job_index_list_forGridssAssemble

	output:
	tuple val(tumor_normal_sample_id), path(assembly_bam_per_index_working_dir) into split_assemble_output_forGridssMergeAssembly

	when:
	params.gridss == "on"

	script:
	assembly_bam = "${tumor_normal_sample_id}.assembly.bam"
	assembly_bam_per_index_working_dir = "${assembly_bam}.gridss.working_${index}"
	"""
	touch gridss.properties
	echo "chunkSize=5000000" >> gridss.properties

	gridss.sh --steps assemble \
	--jvmheap "${task.memory.toGiga()}"g \
	--otherjvmheap "${task.memory.toGiga()}"g \
	--threads "${task.cpus}" \
	--configuration gridss.properties \
	--picardoptions TMP_DIR=. \
	--picardoptions MAX_RECORDS_IN_RAM=4000000 \
	--jobindex "${index}" \
	--jobnodes 2 \
	--assembly "${assembly_bam}" \
	--reference "${reference_genome_fasta_forGridssSetup}" \
	--blacklist "${gatk_bundle_wgs_bed_blacklist_1based_forGridssSetup}" \
	"${normal_bam}" "${tumor_bam}"

	mv "${assembly_bam}.gridss.working" "${assembly_bam_per_index_working_dir}"
	"""
}

// GRIDSS assemble ~ merge all assembly results together
process mergeAssembly_gridss {
	publishDir "${params.output_dir}/somatic/gridss/intermediates", mode: 'symlink', pattern: '*${assembly_bam}'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(assembly_bam_per_index_working_dir), val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(reference_genome_fasta_forGridssSetup), path(reference_genome_fasta_index_forGridssSetup), path(gatk_bundle_wgs_bed_blacklist_1based_forGridssSetup), path("Homo_sapiens_assembly38.fasta.amb"), path("Homo_sapiens_assembly38.fasta.ann"), path("Homo_sapiens_assembly38.fasta.bwt"), path("Homo_sapiens_assembly38.fasta.dict"), path("Homo_sapiens_assembly38.fasta.gridsscache"), path("Homo_sapiens_assembly38.fasta.img"), path("Homo_sapiens_assembly38.fasta.pac"), path("Homo_sapiens_assembly38.fasta.sa"), path(tumor_bam_working_dir), path(normal_bam_working_dir) from split_assemble_output_forGridssMergeAssembly.groupTuple().combine(preprocess_output_forGridssMergeAssembly)

	output:
	tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(reference_genome_fasta_forGridssSetup), path(reference_genome_fasta_index_forGridssSetup), path(gatk_bundle_wgs_bed_blacklist_1based_forGridssSetup), path("Homo_sapiens_assembly38.fasta.amb"), path("Homo_sapiens_assembly38.fasta.ann"), path("Homo_sapiens_assembly38.fasta.bwt"), path("Homo_sapiens_assembly38.fasta.dict"), path("Homo_sapiens_assembly38.fasta.gridsscache"), path("Homo_sapiens_assembly38.fasta.img"), path("Homo_sapiens_assembly38.fasta.pac"), path("Homo_sapiens_assembly38.fasta.sa"), path(tumor_bam_working_dir), path(normal_bam_working_dir), path(assembly_bam), path(assembly_bam_working_dir) into merge_assembly_output_forGridssCalling

	when:
	params.gridss == "on"

	script:
	assembly_bam = "${tumor_normal_sample_id}.assembly.bam"
	assembly_bam_working_dir = "${assembly_bam}.gridss.working"
	"""
	touch gridss.properties
	echo "chunkSize=5000000" >> gridss.properties

	gridss.sh --steps assemble \
	--jvmheap "${task.memory.toGiga() - 4}"g \
	--otherjvmheap "${task.memory.toGiga() - 4}"g \
	--threads "${task.cpus}" \
	--configuration gridss.properties \
	--picardoptions TMP_DIR=. \
	--picardoptions MAX_RECORDS_IN_RAM=4000000 \
	--assembly "${assembly_bam}" \
	--reference "${reference_genome_fasta_forGridssSetup}" \
	--blacklist "${gatk_bundle_wgs_bed_blacklist_1based_forGridssSetup}" \
	"${normal_bam}" "${tumor_bam}"
	"""
}

// GRIDSS call ~ calls structural variants based on alignment-guided positional de Bruijn graph genome
process svAndIndelCalling_gridss {
	publishDir "${params.output_dir}/somatic/gridss/rawVcfs", mode: 'symlink', pattern: '*.{vcf.gz,tbi}'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(reference_genome_fasta_forGridssSetup), path(reference_genome_fasta_index_forGridssSetup), path(gatk_bundle_wgs_bed_blacklist_1based_forGridssSetup), path("Homo_sapiens_assembly38.fasta.amb"), path("Homo_sapiens_assembly38.fasta.ann"), path("Homo_sapiens_assembly38.fasta.bwt"), path("Homo_sapiens_assembly38.fasta.dict"), path("Homo_sapiens_assembly38.fasta.gridsscache"), path("Homo_sapiens_assembly38.fasta.img"), path("Homo_sapiens_assembly38.fasta.pac"), path("Homo_sapiens_assembly38.fasta.sa"), path(tumor_bam_working_dir), path(normal_bam_working_dir), path(assembly_bam), path(assembly_bam_working_dir) from merge_assembly_output_forGridssCalling

	output:
	tuple val(tumor_normal_sample_id), val(tumor_id), val(normal_id), path(raw_gridss_vcf), path(raw_gridss_vcf_index) into raw_vcf_forGridssPostprocessing

	when:
	params.gridss == "on"

	script:
	tumor_id = "${tumor_bam.baseName}".replaceFirst(/\..*$/, "")
	normal_id = "${normal_bam.baseName}".replaceFirst(/\..*$/, "")
	raw_gridss_vcf = "${tumor_normal_sample_id}.gridss.raw.vcf.gz"
	raw_gridss_vcf_index = "${raw_gridss_vcf}.tbi"
	"""
	touch gridss.properties
	echo "chunkSize=5000000" >> gridss.properties

	gridss.sh --steps call \
	--jvmheap "${task.memory.toGiga() - 4}"g \
	--otherjvmheap "${task.memory.toGiga() - 4}"g \
	--threads "${task.cpus}" \
	--configuration gridss.properties \
	--picardoptions TMP_DIR=. \
	--picardoptions MAX_RECORDS_IN_RAM=4000000 \
	--assembly "${assembly_bam}" \
	--reference "${reference_genome_fasta_forGridssSetup}" \
	--blacklist "${gatk_bundle_wgs_bed_blacklist_1based_forGridssSetup}" \
	--output "${tumor_normal_sample_id}" \
	"${normal_bam}" "${tumor_bam}"

	bgzip < "${tumor_normal_sample_id}.gridss.working/${tumor_normal_sample_id}.allocated.vcf" > "${raw_gridss_vcf}"
	tabix "${raw_gridss_vcf}"
	"""
}

// Combine Combine all needed reference FASTA and PoN SV BED/BEDPE files into one channel for GRIDSS setupreference process
reference_genome_fasta_forGridssPostprocessing.combine( reference_genome_fasta_index_forGridssPostprocessing )
	.combine( reference_genome_fasta_dict_forGridssPostprocessing )
	.set{ reference_genome_bundle_forGridssPostprocessing }

reference_genome_bundle_forGridssPostprocessing.combine( pon_single_breakend_bed )
	.combine( pon_breakpoint_bedpe )
	.combine( known_fusion_pairs_bedpe )
	.set{ ref_genome_and_pon_beds_forGridssPostprocessing }

// GRIPSS ~ apply a set of filtering and post processing steps to raw GRIDSS calls
process filteringAndPostprocessesing_gridss {
	publishDir "${params.output_dir}/somatic/gridss", mode: 'copy', pattern: '*.{vcf.gz,tbi}'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), val(tumor_id), val(normal_id), path(raw_gridss_vcf), path(raw_gridss_vcf_index), path(reference_genome_fasta_forGridssPostprocessing), path(reference_genome_fasta_index_forGridssPostprocessing), path(reference_genome_fasta_dict_forGridssPostprocessing), path(pon_single_breakend_bed), path(pon_breakpoint_bedpe), path(known_fusion_pairs_bedpe) from raw_vcf_forGridssPostprocessing.combine(ref_genome_and_pon_beds_forGridssPostprocessing)

	output:
	tuple path(final_gridss_vcf), path(final_gridss_vcf_index)

	when:
	params.gridss == "on"

	script:
	intermediate_filterted_vcf = "${raw_gridss_vcf}".replaceFirst(/\.raw\.vcf\.gz/, ".intermediate.vcf.gz")
	filterted_vcf = "${raw_gridss_vcf}".replaceFirst(/\.raw\.vcf\.gz/, "filtered.vcf.gz")
	final_gridss_vcf = "${raw_gridss_vcf}".replaceFirst(/\.raw\.vcf\.gz/, ".somatic.sv.vcf.gz")
	final_gridss_vcf_index = "${final_gridss_vcf}.tbi"
	"""
	java -Xmx${task.memory.toGiga()}G -cp /opt/gripss/gripss.jar com.hartwig.hmftools.gripss.GripssApplicationKt \
	-tumor "${tumor_id}" \
    -reference "${normal_id}" \
    -ref_genome "${reference_genome_fasta_forGridssPostprocessing}" \
    -breakend_pon "${pon_single_breakend_bed}" \
    -breakpoint_pon "${pon_breakpoint_bedpe}" \
    -breakpoint_hotspot "${known_fusion_pairs_bedpe}" \
    -input_vcf "${raw_gridss_vcf}" \
    -output_vcf "${intermediate_filterted_vcf}"

    java -Xmx${task.memory.toGiga()}G -cp /opt/gripss/gripss.jar com.hartwig.hmftools.gripss.GripssHardFilterApplicationKt \
    -input_vcf  "${intermediate_filterted_vcf}" \
    -output_vcf "${filterted_vcf}"

    zgrep -E '^#|PASS' "${filterted_vcf}" \
    | \
    bgzip > "${final_gridss_vcf}"
    tabix "${final_gridss_vcf}"
	"""
}

// END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \\




/*


// ~~~~~~~~~ CONSENSUS SNV/INDEL ~~~~~~~~~~~ \\
// START

// MergeVCF ~ merge VCF files by calls, labelling calls by the callers that made them to generate consensus
process mergeAndGenerateConsensusSnvCalls_mergevcf {
	publishDir "${params.output_dir}/somatic/consensus/intermediates", mode: 'copy', pattern: '*.{vcf.gz,tbi}'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(final_varscan_snv_vcf), path(final_varscan_snv_vcf_index), path(final_mutect_snv_vcf), path(final_mutect_snv_vcf_index), path(final_caveman_snv_vcf), path(final_caveman_snv_vcf_index), path(final_strelka_snv_vcf), path(final_strelka_snv_vcf_index) from final_varscan_snv_vcf_forConsensus.join(final_mutect_snv_vcf_forConsensus).join(final_caveman_snv_vcf_forConsensus).join(final_strelka_snv_vcf_forConsensus)

	output:
	tuple val(tumor_normal_sample_id), path(consensus_somatic_snv_nosamples_badheader_noformat_vcf), path(snv_vcf_base_header) into consensus_snv_vcf_forConsensusSnvMpileup

	when:
	params.varscan == "on" && params.mutect == "on" && params.strelka == "on" && params.caveman == "on"

	script:
	consensus_somatic_snv_nosamples_badheader_noformat_vcf = "${tumor_normal_sample_id}.consensus.somatic.snv.nosamples.badheader.noformat.vcf"
	snv_vcf_base_header = "snv_vcf_base_header.txt"
	"""
	mergevcf \
	--labels varscan,mutect,strelka,caveman \
	--ncallers \
	--mincallers 2 \
	"${final_varscan_snv_vcf}" \
	"${final_mutect_snv_vcf}" \
	"${final_caveman_snv_vcf}" \
	"${final_strelka_snv_vcf}" \
	| \
	awk '\$1 ~ /^#/ {print \$0;next} {print \$0 | "sort -k1,1V -k2,2n"}' > "${consensus_somatic_snv_nosamples_badheader_noformat_vcf}"

	touch "${snv_vcf_base_header}"
	grep '##fileformat=' "${consensus_somatic_snv_nosamples_badheader_noformat_vcf}" >> "${snv_vcf_base_header}"
	zgrep '##fileDate=' "${final_strelka_snv_vcf}" >> "${snv_vcf_base_header}"
	echo '##source=varscan,mutect,strelka,caveman' >> "${snv_vcf_base_header}"
	zgrep '##normal_sample=' "${final_mutect_snv_vcf}" >> "${snv_vcf_base_header}"
	zgrep '##tumor_sample=' "${final_mutect_snv_vcf}" >> "${snv_vcf_base_header}"
	zgrep '##reference=' "${final_caveman_snv_vcf}" >> "${snv_vcf_base_header}"
	zgrep '##contig=<ID=chr' "${final_mutect_snv_vcf}" >> "${snv_vcf_base_header}"
	"""
}

// MergeVCF ~ merge VCF files by calls, labelling calls by the callers that made them to generate consensus
process mergeAndGenerateConsensusIndelCalls_mergevcf {
	publishDir "${params.output_dir}/somatic/consensus/intermediates", mode: 'copy', pattern: '*.{vcf.gz,tbi}'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(final_varscan_indel_vcf), path(final_varscan_indel_vcf_index), path(final_mutect_indel_vcf), path(final_mutect_indel_vcf_index), path(final_strelka_indel_vcf), path(final_strelka_indel_vcf_index), path(final_svaba_indel_vcf), path(final_svaba_indel_vcf_index) from final_varscan_indel_vcf_forConsensus.join(final_mutect_indel_vcf_forConsensus).join(final_strelka_indel_vcf_forConsensus).join(final_svaba_indel_vcf_forConsensus)

	output:
	tuple val(tumor_normal_sample_id), path(consensus_somatic_indel_nosamples_badheader_noformat_vcf), path(indel_vcf_base_header) into consensus_indel_vcf_forConsensusIndelMpileup

	when:
	params.varscan == "on" && params.mutect == "on" && params.strelka == "on" && params.svaba == "on"

	script:
	consensus_somatic_indel_nosamples_badheader_noformat_vcf = "${tumor_normal_sample_id}.consensus.somatic.indel.nosamples.badheader.noformat.vcf"
	indel_vcf_base_header = "indel_vcf_base_header.txt"
	"""
	mergevcf \
	--labels varscan,mutect,strelka,svaba \
	--ncallers \
	--mincallers 2 \
	"${final_varscan_indel_vcf}" \
	"${final_mutect_indel_vcf}" \
	"${final_strelka_indel_vcf}" \
	"${final_svaba_indel_vcf}" \
	| \
	awk '\$1 ~ /^#/ {print \$0;next} {print \$0 | "sort -k1,1V -k2,2n"}' > "${consensus_somatic_indel_nosamples_badheader_noformat_vcf}"

	touch "${indel_vcf_base_header}"
	grep '##fileformat=' "${consensus_somatic_indel_nosamples_badheader_noformat_vcf}" >> "${indel_vcf_base_header}"
	zgrep '##fileDate=' "${final_svaba_indel_vcf}" >> "${indel_vcf_base_header}"
	echo '##source=varscan,mutect,strelka,svaba' >> "${indel_vcf_base_header}"
	zgrep '##normal_sample=' "${final_mutect_indel_vcf}" >> "${indel_vcf_base_header}"
	zgrep '##tumor_sample=' "${final_mutect_indel_vcf}" >> "${indel_vcf_base_header}"
	zgrep '##reference=' "${final_svaba_indel_vcf}" >> "${indel_vcf_base_header}"
	zgrep '##contig=<ID=chr' "${final_mutect_indel_vcf}" >> "${indel_vcf_base_header}"
	"""
}

// END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \\


// ~~~~~~ COMPLETE CONSENSUS SNV VCF ~~~~~~~ \\
// START

// Combine all needed reference FASTA files into one channel for use in BCFtools mpileup
reference_genome_fasta_forConsensusSnvMpileup.combine( reference_genome_fasta_index_forConsensusSnvMpileup )
	.combine( reference_genome_fasta_dict_forConsensusSnvMpileup )
	.set{ reference_genome_bundle_forConsensusSnvMpileup }

// BCFtools mpileup ~ generate mpileup metrics for consensus SNV variant calls
process consensusSnvMpileup_bcftools {
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(consensus_somatic_snv_nosamples_badheader_noformat_vcf), path(snv_vcf_base_header), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(reference_genome_fasta_forConsensusSnvMpileup), path(reference_genome_fasta_index_forConsensusSnvMpileup), path(reference_genome_fasta_dict_forConsensusSnvMpileup) from consensus_snv_vcf_forConsensusSnvMpileup.join(bams_forConsensusSnvMpileup).combine(reference_genome_bundle_forConsensusSnvMpileup)

	output:
	tuple val(tumor_normal_sample_id), val(tumor_id), val(normal_id), path(mpileup_supported_consensus_somatic_snv_nosamples_noformat_vcf) into consensus_snv_vcf_forAddSamples
	tuple val(tumor_normal_sample_id), path(snv_mpileup_info_dp_metrics), path(snv_mpileup_normal_format_metrics), path(snv_mpileup_normal_format_metrics_index), path(snv_mpileup_tumor_format_metrics), path(snv_mpileup_tumor_format_metrics_index) into consensus_snv_mpileup_metrics_forAddFormat

	when:
	params.varscan == "on" && params.mutect == "on" && params.strelka == "on" && params.caveman == "on"

	script:
	tumor_id = "${tumor_bam.baseName}".replaceFirst(/\..*$/, "")
	normal_id = "${normal_bam.baseName}".replaceFirst(/\..*$/, "")
	full_snv_vcf_header = "full_snv_vcf_header.txt"
	mpileup_supported_consensus_somatic_snv_nosamples_noformat_vcf = "${tumor_normal_sample_id}.ms.consensus.somatic.snv.nosamples.noformat.vcf"
	snv_mpileup_info_dp_metrics = "${tumor_normal_sample_id}.snv.mpileup.info.dp.metrics.txt"
	snv_mpileup_normal_format_metrics = "${tumor_normal_sample_id}.snv.mpileup.normal.format.metrics.txt.gz"
	snv_mpileup_normal_format_metrics_index = "${snv_mpileup_normal_format_metrics}.tbi"
	snv_mpileup_tumor_format_metrics = "${tumor_normal_sample_id}.snv.mpileup.tumor.format.metrics.txt.gz"
	snv_mpileup_tumor_format_metrics_index = "${snv_mpileup_tumor_format_metrics}.tbi"
	"""
	bcftools mpileup \
	--no-BAQ \
	--min-MQ 35 \
	--min-BQ 30 \
	--max-depth 5000 \
	--threads ${task.cpus} \
	--fasta-ref "${reference_genome_fasta_forConsensusSnvMpileup}" \
	--regions-file "${consensus_somatic_snv_nosamples_badheader_noformat_vcf}" \
	--samples ${normal_id},${tumor_id} \
	--annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP \
	"${normal_bam}" "${tumor_bam}" \
	| \
	bcftools norm \
	--threads ${task.cpus} \
	--multiallelics - \
	--fasta-ref "${reference_genome_fasta_forConsensusSnvMpileup}" \
	- \
	| \
	bcftools filter \
	--exclude 'FORMAT/DP == 0' \
	- \
	| \
	grep -Ev '<.>|INDEL;' \
	| \
	bgzip > "${tumor_normal_sample_id}.consensus.somatic.snv.mpileup.vcf.gz"
	tabix "${tumor_normal_sample_id}.consensus.somatic.snv.mpileup.vcf.gz"

	bgzip < "${consensus_somatic_snv_nosamples_badheader_noformat_vcf}" > "${consensus_somatic_snv_nosamples_badheader_noformat_vcf}.gz"
	tabix "${consensus_somatic_snv_nosamples_badheader_noformat_vcf}.gz"

	touch "${full_snv_vcf_header}"
	cat "${snv_vcf_base_header}" >> "${full_snv_vcf_header}"
	zgrep '##FILTER=<ID=LOWSUPPORT' "${consensus_somatic_snv_nosamples_badheader_noformat_vcf}.gz" >> "${full_snv_vcf_header}"
	zgrep '##INFO=' "${consensus_somatic_snv_nosamples_badheader_noformat_vcf}.gz" >> "${full_snv_vcf_header}"
	echo -e '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO' >> "${full_snv_vcf_header}"

	bcftools isec \
	--nfiles =2 \
	--write 1 \
	"${consensus_somatic_snv_nosamples_badheader_noformat_vcf}.gz" \
	"${tumor_normal_sample_id}.consensus.somatic.snv.mpileup.vcf.gz" \
	| \
	bcftools reheader \
	--header "${full_snv_vcf_header}" \
	--output "${mpileup_supported_consensus_somatic_snv_nosamples_noformat_vcf}"

	bcftools query \
	--format '%CHROM\t%POS\t%REF\t%ALT\t%INFO/DP\n' \
	--output "${snv_mpileup_info_dp_metrics}" \
	"${tumor_normal_sample_id}.consensus.somatic.snv.mpileup.vcf.gz"

	bcftools query \
	--format '%CHROM\t%POS\t%REF\t%ALT\t[%DP]\t[%AD]\t[%ADF]\t[%ADR]\n' \
	--samples "${normal_id}" \
	"${tumor_normal_sample_id}.consensus.somatic.snv.mpileup.vcf.gz" \
	| \
	bgzip > "${snv_mpileup_normal_format_metrics}"
	tabix -s1 -b2 -e2 "${snv_mpileup_normal_format_metrics}"

	bcftools query \
	--format '%CHROM\t%POS\t%REF\t%ALT\t[%DP]\t[%AD]\t[%ADF]\t[%ADR]\n' \
	--samples "${tumor_id}" \
	"${tumor_normal_sample_id}.consensus.somatic.snv.mpileup.vcf.gz" \
	| \
	bgzip > "${snv_mpileup_tumor_format_metrics}"
	tabix -s1 -b2 -e2 "${snv_mpileup_tumor_format_metrics}"
	"""
}

// VAtools vcf-genotype-annotator ~ add samples to VCF and fill in placeholder genotype FORMAT field
process addSamplesToConsensusSnvVcf_vatools {
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), val(tumor_id), val(normal_id), path(mpileup_supported_consensus_somatic_snv_nosamples_noformat_vcf) from consensus_snv_vcf_forAddSamples

	output:
	tuple val(tumor_normal_sample_id), val(tumor_id), val(normal_id), path(mpileup_supported_consensus_somatic_snv_noformat_vcf) into consensus_snv_vcf_forAddFormat

	when:
	params.varscan == "on" && params.mutect == "on" && params.strelka == "on" && params.caveman == "on"

	script:
	mpileup_supported_consensus_somatic_snv_noformat_vcf = "${tumor_normal_sample_id}.ms.consensus.somatic.snv.noformat.vcf"
	"""
	vcf-genotype-annotator \
	--output-vcf "${tumor_normal_sample_id}.ms.consensus.somatic.snv.halfsamples.noformat.vcf" \
	"${mpileup_supported_consensus_somatic_snv_nosamples_noformat_vcf}" \
	"${normal_id}" \
	.

	vcf-genotype-annotator \
	--output-vcf "${mpileup_supported_consensus_somatic_snv_noformat_vcf}" \
	"${tumor_normal_sample_id}.ms.consensus.somatic.snv.halfsamples.noformat.vcf" \
	"${tumor_id}" \
	.
	"""
}

// BCFtools annotate ~ modify VCF INFO/FORMAT columns to include better information and final filtering
process annotateConsensusSnvVcfFormatColumnAndFilter_bcftools {
	publishDir "${params.output_dir}/somatic/consensus", mode: 'copy', pattern: '*.{vcf.gz,tbi}'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), val(tumor_id), val(normal_id), path(mpileup_supported_consensus_somatic_snv_noformat_vcf), path(snv_mpileup_info_dp_metrics), path(snv_mpileup_normal_format_metrics), path(snv_mpileup_normal_format_metrics_index), path(snv_mpileup_tumor_format_metrics), path(snv_mpileup_tumor_format_metrics_index) from consensus_snv_vcf_forAddFormat.join(consensus_snv_mpileup_metrics_forAddFormat)

	output:
	tuple path(hq_snv_consensus_vcf), path(hq_snv_consensus_vcf_index) into high_quality_consensus_snv_forAnnotation

	when:
	params.varscan == "on" && params.mutect == "on" && params.strelka == "on" && params.caveman == "on"

	script:
	hq_snv_consensus_vcf_info_header = "hq_snv_consensus_vcf_info_header.txt"
	hq_snv_consensus_vcf_format_headers = "hq_snv_consensus_vcf_format_header.txt"
	hq_snv_consensus_vcf = "${tumor_normal_sample_id}.hq.consensus.somatic.snv.vcf.gz"
	hq_snv_consensus_vcf_index = "${hq_snv_consensus_vcf}.tbi"
	"""
	cat "${snv_mpileup_info_dp_metrics}" \
	| \
	paste - <(zcat "${snv_mpileup_tumor_format_metrics}" | cut -f 6 | awk '{split(\$0,x,","); print x[2]}') > "${tumor_normal_sample_id}.snv.mpileup.info.dp.ac.metrics.txt"

	cat "${tumor_normal_sample_id}.snv.mpileup.info.dp.ac.metrics.txt" \
	| \
	paste - <(zcat "${snv_mpileup_tumor_format_metrics}" | cut -f 5) \
	| \
	awk 'BEGIN {OFS="\t"} {print \$1,\$2,\$3,\$4,\$5,\$6,\$6/\$7}' \
	| \
	bgzip > "${tumor_normal_sample_id}.snv.mpileup.info.metrics.txt.gz"
	tabix -s1 -b2 -e2 "${tumor_normal_sample_id}.snv.mpileup.info.metrics.txt.gz"

	touch "${hq_snv_consensus_vcf_info_header}"
	echo '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth across samples (normal sample DP + tumor sample DP)">' >> "${hq_snv_consensus_vcf_info_header}"
	echo '##INFO=<ID=AC,Number=1,Type=Integer,Description="Count of ALT allele reads in tumor sample">' >> "${hq_snv_consensus_vcf_info_header}"
	echo '##INFO=<ID=VAF,Number=1,Type=Float,Description="Variant allele frequency, expressed as fraction of ALT allele reads in total read depth in tumor sample (tumor sample ALT AC / tumor sample DP)">' >> "${hq_snv_consensus_vcf_info_header}"

	bcftools annotate \
	--output-type z \
	--annotations "${tumor_normal_sample_id}.snv.mpileup.info.metrics.txt.gz" \
	--header-lines "${hq_snv_consensus_vcf_info_header}" \
	--columns CHROM,POS,REF,ALT,INFO/DP,INFO/AC,INFO/VAF \
	--output "${tumor_normal_sample_id}.ms.consensus.somatic.snv.info.noformat.vcf.gz" \
	"${mpileup_supported_consensus_somatic_snv_noformat_vcf}"

	touch "${hq_snv_consensus_vcf_format_headers}"
	echo '##FORMAT=<ID=DPS,Number=1,Type=Integer,Description="Total read depth in sample">' >> "${hq_snv_consensus_vcf_format_headers}"
	echo '##FORMAT=<ID=ACS,Number=R,Type=Integer,Description="Count of REF,ALT allele reads in sample">' >> "${hq_snv_consensus_vcf_format_headers}"
	echo '##FORMAT=<ID=ACFS,Number=R,Type=Integer,Description="Count of REF,ALT allele reads on forward(+) strand in sample">' >> "${hq_snv_consensus_vcf_format_headers}"
	echo '##FORMAT=<ID=ACRS,Number=R,Type=Integer,Description="Count of REF,ALT allele reads on reverse(-) strand in sample">' >> "${hq_snv_consensus_vcf_format_headers}"

	bcftools annotate \
	--output-type z \
	--samples "${normal_id}" \
	--annotations "${snv_mpileup_normal_format_metrics}" \
	--header-lines "${hq_snv_consensus_vcf_format_headers}" \
	--columns CHROM,POS,REF,ALT,FORMAT/DPS,FORMAT/ACS,FORMAT/ACFS,FORMAT/ACRS \
	--output "${tumor_normal_sample_id}.ms.consensus.somatic.snv.info.halfformat.vcf.gz" \
	"${tumor_normal_sample_id}.ms.consensus.somatic.snv.info.noformat.vcf.gz"

	bcftools annotate \
	--output-type z \
	--samples "${tumor_id}" \
	--annotations "${snv_mpileup_tumor_format_metrics}" \
	--header-lines "${hq_snv_consensus_vcf_format_headers}" \
	--columns CHROM,POS,REF,ALT,FORMAT/DPS,FORMAT/ACS,FORMAT/ACFS,FORMAT/ACRS \
	--remove FORMAT/GT \
	--output "${tumor_normal_sample_id}.ms.consensus.somatic.snv.info.format.vcf.gz" \
	"${tumor_normal_sample_id}.ms.consensus.somatic.snv.info.halfformat.vcf.gz"

	bcftools filter \
	--output-type v \
	--exclude 'INFO/AC<2 | INFO/VAF<0.01' \
	"${tumor_normal_sample_id}.ms.consensus.somatic.snv.info.format.vcf.gz" \
	| \
	grep -v '##bcftools_annotate' \
	| \
	bgzip > "${hq_snv_consensus_vcf}"
	tabix "${hq_snv_consensus_vcf}"
	"""
}

// END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \\


// ~~~~~ COMPLETE CONSENSUS INDEL VCF ~~~~~~ \\
// START

// Combine all needed reference FASTA files into one channel for use in BCFtools mpileup
reference_genome_fasta_forConsensusIndelMpileup.combine( reference_genome_fasta_index_forConsensusIndelMpileup )
	.combine( reference_genome_fasta_dict_forConsensusIndelMpileup )
	.set{ reference_genome_bundle_forConsensusIndelMpileup }

// BCFtools mpileup ~ generate mpileup metrics for consensus indel variant calls
process consensusIndelMpileup_bcftools {
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(consensus_somatic_indel_nosamples_badheader_noformat_vcf), path(indel_vcf_base_header), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(reference_genome_fasta_forConsensusIndelMpileup), path(reference_genome_fasta_index_forConsensusIndelMpileup), path(reference_genome_fasta_dict_forConsensusIndelMpileup) from consensus_indel_vcf_forConsensusIndelMpileup.join(bams_forConsensusIndelMpileup).combine(reference_genome_bundle_forConsensusIndelMpileup)

	output:
	tuple val(tumor_normal_sample_id), val(tumor_id), val(normal_id), path(mpileup_supported_consensus_somatic_indel_nosamples_noformat_vcf) into consensus_indel_vcf_forAddSamples
	tuple val(tumor_normal_sample_id), path(indel_mpileup_info_dp_metrics), path(indel_mpileup_normal_format_metrics), path(indel_mpileup_normal_format_metrics_index), path(indel_mpileup_tumor_format_metrics), path(indel_mpileup_tumor_format_metrics_index) into consensus_indel_mpileup_metrics_forAddFormat

	when:
	params.varscan == "on" && params.mutect == "on" && params.strelka == "on" && params.svaba == "on"

	script:
	tumor_id = "${tumor_bam.baseName}".replaceFirst(/\..*$/, "")
	normal_id = "${normal_bam.baseName}".replaceFirst(/\..*$/, "")
	full_indel_vcf_header = "full_indel_vcf_header.txt"
	mpileup_supported_consensus_somatic_indel_nosamples_noformat_vcf = "${tumor_normal_sample_id}.ms.consensus.somatic.indel.nosamples.noformat.vcf"
	indel_mpileup_info_dp_metrics = "${tumor_normal_sample_id}.indel.mpileup.info.dp.metrics.txt"
	indel_mpileup_normal_format_metrics = "${tumor_normal_sample_id}.indel.mpileup.normal.format.metrics.txt.gz"
	indel_mpileup_normal_format_metrics_index = "${indel_mpileup_normal_format_metrics}.tbi"
	indel_mpileup_tumor_format_metrics = "${tumor_normal_sample_id}.indel.mpileup.tumor.format.metrics.txt.gz"
	indel_mpileup_tumor_format_metrics_index = "${indel_mpileup_tumor_format_metrics}.tbi"
	"""
	bcftools mpileup \
	--no-BAQ \
	--min-MQ 35 \
	--min-BQ 30 \
	--max-depth 5000 \
	--threads ${task.cpus} \
	--fasta-ref "${reference_genome_fasta_forConsensusIndelMpileup}" \
	--regions-file "${consensus_somatic_indel_nosamples_badheader_noformat_vcf}" \
	--samples ${normal_id},${tumor_id} \
	--annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP \
	"${normal_bam}" "${tumor_bam}" \
	| \
	bcftools norm \
	--threads ${task.cpus} \
	--multiallelics - \
	--fasta-ref "${reference_genome_fasta_forConsensusIndelMpileup}" \
	- \
	| \
	bcftools filter \
	--exclude 'FORMAT/DP == 0' \
	- \
	| \
	grep -E '^#|INDEL;' \
	| \
	bgzip > "${tumor_normal_sample_id}.consensus.somatic.indel.mpileup.vcf.gz"
	tabix "${tumor_normal_sample_id}.consensus.somatic.indel.mpileup.vcf.gz"

	bgzip < "${consensus_somatic_indel_nosamples_badheader_noformat_vcf}" > "${consensus_somatic_indel_nosamples_badheader_noformat_vcf}.gz"
	tabix "${consensus_somatic_indel_nosamples_badheader_noformat_vcf}.gz"

	touch "${full_indel_vcf_header}"
	cat "${indel_vcf_base_header}" >> "${full_indel_vcf_header}"
	zgrep '##FILTER=<ID=LOWSUPPORT' "${consensus_somatic_indel_nosamples_badheader_noformat_vcf}.gz" >> "${full_indel_vcf_header}"
	zgrep '##INFO=' "${consensus_somatic_indel_nosamples_badheader_noformat_vcf}.gz" >> "${full_indel_vcf_header}"
	echo -e '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO' >> "${full_indel_vcf_header}"

	bcftools isec \
	--nfiles =2 \
	--write 1 \
	"${consensus_somatic_indel_nosamples_badheader_noformat_vcf}.gz" \
	"${tumor_normal_sample_id}.consensus.somatic.indel.mpileup.vcf.gz" \
	| \
	bcftools reheader \
	--header "${full_indel_vcf_header}" \
	--output "${mpileup_supported_consensus_somatic_indel_nosamples_noformat_vcf}"

	bcftools query \
	--format '%CHROM\t%POS\t%REF\t%ALT\t%INFO/DP\n' \
	--output "${indel_mpileup_info_dp_metrics}" \
	"${tumor_normal_sample_id}.consensus.somatic.indel.mpileup.vcf.gz"

	bcftools query \
	--format '%CHROM\t%POS\t%REF\t%ALT\t[%DP]\t[%AD]\t[%ADF]\t[%ADR]\n' \
	--samples "${normal_id}" \
	"${tumor_normal_sample_id}.consensus.somatic.indel.mpileup.vcf.gz" \
	| \
	bgzip > "${indel_mpileup_normal_format_metrics}"
	tabix -s1 -b2 -e2 "${indel_mpileup_normal_format_metrics}"

	bcftools query \
	--format '%CHROM\t%POS\t%REF\t%ALT\t[%DP]\t[%AD]\t[%ADF]\t[%ADR]\n' \
	--samples "${tumor_id}" \
	"${tumor_normal_sample_id}.consensus.somatic.indel.mpileup.vcf.gz" \
	| \
	bgzip > "${indel_mpileup_tumor_format_metrics}"
	tabix -s1 -b2 -e2 "${indel_mpileup_tumor_format_metrics}"
	"""
}

// VAtools vcf-genotype-annotator ~ add samples to VCF and fill in placeholder genotype FORMAT field
process addSamplesToConsensusIndelVcf_vatools {
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), val(tumor_id), val(normal_id), path(mpileup_supported_consensus_somatic_indel_nosamples_noformat_vcf) from consensus_indel_vcf_forAddSamples

	output:
	tuple val(tumor_normal_sample_id), val(tumor_id), val(normal_id), path(mpileup_supported_consensus_somatic_indel_noformat_vcf) into consensus_indel_vcf_forAddFormat

	when:
	params.varscan == "on" && params.mutect == "on" && params.strelka == "on" && params.svaba == "on"

	script:
	mpileup_supported_consensus_somatic_indel_noformat_vcf = "${tumor_normal_sample_id}.ms.consensus.somatic.indel.noformat.vcf"
	"""
	vcf-genotype-annotator \
	--output-vcf "${tumor_normal_sample_id}.ms.consensus.somatic.indel.halfsamples.noformat.vcf" \
	"${mpileup_supported_consensus_somatic_indel_nosamples_noformat_vcf}" \
	"${normal_id}" \
	.

	vcf-genotype-annotator \
	--output-vcf "${mpileup_supported_consensus_somatic_indel_noformat_vcf}" \
	"${tumor_normal_sample_id}.ms.consensus.somatic.indel.halfsamples.noformat.vcf" \
	"${tumor_id}" \
	.
	"""
}

// BCFtools annotate ~ modify VCF INFO/FORMAT columns to include better information and final filtering
process annotateConsensusIndelVcfFormatColumnAndFilter_bcftools {
	publishDir "${params.output_dir}/somatic/consensus", mode: 'copy', pattern: '*.{vcf.gz,tbi}'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), val(tumor_id), val(normal_id), path(mpileup_supported_consensus_somatic_indel_noformat_vcf), path(indel_mpileup_info_dp_metrics), path(indel_mpileup_normal_format_metrics), path(indel_mpileup_normal_format_metrics_index), path(indel_mpileup_tumor_format_metrics), path(indel_mpileup_tumor_format_metrics_index) from consensus_indel_vcf_forAddFormat.join(consensus_indel_mpileup_metrics_forAddFormat)

	output:
	tuple path(hq_indel_consensus_vcf), path(hq_indel_consensus_vcf_index) into high_quality_consensus_indel_forAnnotation

	when:
	params.varscan == "on" && params.mutect == "on" && params.strelka == "on" && params.svaba == "on"

	script:
	hq_indel_consensus_vcf_info_header = "hq_indel_consensus_vcf_info_header.txt"
	hq_indel_consensus_vcf_format_headers = "hq_indel_consensus_vcf_format_header.txt"
	hq_indel_consensus_vcf = "${tumor_normal_sample_id}.hq.consensus.somatic.indel.vcf.gz"
	hq_indel_consensus_vcf_index = "${hq_indel_consensus_vcf}.tbi"
	"""
	cat "${indel_mpileup_info_dp_metrics}" \
	| \
	paste - <(zcat "${indel_mpileup_tumor_format_metrics}" | cut -f 6 | awk '{split(\$0,x,","); print x[2]}') > "${tumor_normal_sample_id}.indel.mpileup.info.dp.ac.metrics.txt"

	cat "${tumor_normal_sample_id}.indel.mpileup.info.dp.ac.metrics.txt" \
	| \
	paste - <(zcat "${indel_mpileup_tumor_format_metrics}" | cut -f 5) \
	| \
	awk 'BEGIN {OFS="\t"} {print \$1,\$2,\$3,\$4,\$5,\$6,\$6/\$7}' \
	| \
	bgzip > "${tumor_normal_sample_id}.indel.mpileup.info.metrics.txt.gz"
	tabix -s1 -b2 -e2 "${tumor_normal_sample_id}.indel.mpileup.info.metrics.txt.gz"

	touch "${hq_indel_consensus_vcf_info_header}"
	echo '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth across samples (normal sample DP + tumor sample DP)">' >> "${hq_indel_consensus_vcf_info_header}"
	echo '##INFO=<ID=AC,Number=1,Type=Integer,Description="Count of ALT allele reads in tumor sample">' >> "${hq_indel_consensus_vcf_info_header}"
	echo '##INFO=<ID=VAF,Number=1,Type=Float,Description="Variant allele frequency, expressed as fraction of ALT allele reads in total read depth in tumor sample (tumor sample ALT AC / tumor sample DP)">' >> "${hq_indel_consensus_vcf_info_header}"

	bcftools annotate \
	--output-type z \
	--annotations "${tumor_normal_sample_id}.indel.mpileup.info.metrics.txt.gz" \
	--header-lines "${hq_indel_consensus_vcf_info_header}" \
	--columns CHROM,POS,REF,ALT,INFO/DP,INFO/AC,INFO/VAF \
	--output "${tumor_normal_sample_id}.ms.consensus.somatic.indel.info.noformat.vcf.gz" \
	"${mpileup_supported_consensus_somatic_indel_noformat_vcf}"

	touch "${hq_indel_consensus_vcf_format_headers}"
	echo '##FORMAT=<ID=DPS,Number=1,Type=Integer,Description="Total read depth in sample">' >> "${hq_indel_consensus_vcf_format_headers}"
	echo '##FORMAT=<ID=ACS,Number=R,Type=Integer,Description="Count of REF,ALT allele reads in sample">' >> "${hq_indel_consensus_vcf_format_headers}"
	echo '##FORMAT=<ID=ACFS,Number=R,Type=Integer,Description="Count of REF,ALT allele reads on forward(+) strand in sample">' >> "${hq_indel_consensus_vcf_format_headers}"
	echo '##FORMAT=<ID=ACRS,Number=R,Type=Integer,Description="Count of REF,ALT allele reads on reverse(-) strand in sample">' >> "${hq_indel_consensus_vcf_format_headers}"

	bcftools annotate \
	--output-type z \
	--samples "${normal_id}" \
	--annotations "${indel_mpileup_normal_format_metrics}" \
	--header-lines "${hq_indel_consensus_vcf_format_headers}" \
	--columns CHROM,POS,REF,ALT,FORMAT/DPS,FORMAT/ACS,FORMAT/ACFS,FORMAT/ACRS \
	--output "${tumor_normal_sample_id}.ms.consensus.somatic.indel.info.halfformat.vcf.gz" \
	"${tumor_normal_sample_id}.ms.consensus.somatic.indel.info.noformat.vcf.gz"

	bcftools annotate \
	--output-type z \
	--samples "${tumor_id}" \
	--annotations "${indel_mpileup_tumor_format_metrics}" \
	--header-lines "${hq_indel_consensus_vcf_format_headers}" \
	--columns CHROM,POS,REF,ALT,FORMAT/DPS,FORMAT/ACS,FORMAT/ACFS,FORMAT/ACRS \
	--remove FORMAT/GT \
	--output "${tumor_normal_sample_id}.ms.consensus.somatic.indel.info.format.vcf.gz" \
	"${tumor_normal_sample_id}.ms.consensus.somatic.indel.info.halfformat.vcf.gz"

	bcftools filter \
	--output-type v \
	--exclude 'INFO/AC<2 | INFO/VAF<0.01' \
	"${tumor_normal_sample_id}.ms.consensus.somatic.indel.info.format.vcf.gz" \
	| \
	grep -v '##bcftools_annotate' \
	| \
	bgzip > "${hq_indel_consensus_vcf}"
	tabix "${hq_indel_consensus_vcf}"
	"""
}

// END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \\



// ~~~~~~~~~~~~~~~ SURVIVOR ~~~~~~~~~~~~~~~~ \\
// START

// BCFtools reheader / view ~ rename the samples within SvABA VCF and then remove normal sample from both SvABA, Manta, and DELLY VCFs
process prepSvVcfsForSurvivor_bcftools {
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), val(tumor_id), path(final_manta_somatic_sv_vcf), path(final_manta_somatic_sv_vcf_index), path(final_svaba_somatic_sv_vcf), path(final_svaba_somatic_sv_vcf_index), path(sample_renaming_file), path(final_delly_somatic_sv_vcf), path(final_delly_somatic_sv_vcf_index) from manta_sv_vcf_forSurvivorPrep.join(svaba_sv_vcf_forSurvivorPrep, by: [0,1]).join(delly_sv_vcf_forSurvivorPrep, by: [0,1])

	output:
	tuple val(tumor_normal_sample_id), path(manta_tumor_sample_sv_vcf), path(svaba_tumor_sample_sv_vcf), path(delly_tumor_sample_sv_vcf) into sv_vcfs_forSurvivor

	when:
	params.manta == "on" && params.svaba == "on" && params.delly == "on"

	script:
	manta_tumor_sample_sv_vcf = "${tumor_normal_sample_id}.manta.somatic.sv.tumorsample.vcf"
	svaba_tumor_sample_sv_vcf = "${tumor_normal_sample_id}.svaba.somatic.sv.tumorsample.vcf"
	delly_tumor_sample_sv_vcf = "${tumor_normal_sample_id}.delly.somatic.sv.tumorsample.vcf"
	"""
	bcftools view \
	--output-type v \
	--samples "${tumor_id}" \
	--output-file "${manta_tumor_sample_sv_vcf}" \
	"${final_manta_somatic_sv_vcf}"

	bcftools reheader \
	--samples "${sample_renaming_file}" \
	"${final_svaba_somatic_sv_vcf}" \
	| \
	bcftools view \
	--output-type v \
	--samples "${tumor_id}" \
	--output-file "${svaba_tumor_sample_sv_vcf}"

	bcftools view \
	--output-type v \
	--samples "${tumor_id}" \
	--output-file "${delly_tumor_sample_sv_vcf}" \
	"${final_delly_somatic_sv_vcf}"
	"""
}

// SURVIVOR ~ merge SV VCF files to generate a consensus
process mergeAndGenerateConsensusSvCalls_survivor {
	publishDir "${params.output_dir}/somatic/consensus/intermediates", mode: 'copy', pattern: '*.{vcf}'
	tag	"${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(manta_tumor_sample_sv_vcf), path(svaba_tumor_sample_sv_vcf), path(delly_tumor_sample_sv_vcf) from sv_vcfs_forSurvivor

	output:
	tuple val(tumor_normal_sample_id), path(consensus_somatic_sv_vcf)

	when:
	params.manta == "on" && params.svaba == "on" && params.delly == "on"

	script:
	consensus_somatic_sv_vcf = "${tumor_normal_sample_id}.consensus.somatic.sv.vcf"
	"""
	ls *.vcf > input_vcf_list.txt

	SURVIVOR merge \
	input_vcf_list.txt \
	300 \
	2 \
	0 \
	1 \
	0 \
	51 \
	"${consensus_somatic_sv_vcf}"
	"""
}

// END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \\




*/


// ~~~~~~~~~~~~ VEP ANNOTATION ~~~~~~~~~~~~~ \\
// START

// VEP ~ download the reference files used for VEP annotation, if needed
process downloadVepAnnotationReferences_vep {
	publishDir "references/hg38", mode: 'copy'

	output:
	path cached_ref_dir_vep into vep_ref_dir_fromProcess

	when:
	params.vep_ref_cached == "no"

	script:
	cached_ref_dir_vep = "homo_sapiens_vep_101_GRCh38"
	"""
	curl -O ftp://ftp.ensembl.org/pub/release-101/variation/indexed_vep_cache/homo_sapiens_vep_101_GRCh38.tar.gz && \
	mkdir -p "${cached_ref_dir_vep}" && \
	mv homo_sapiens_vep_101_GRCh38.tar.gz "${cached_ref_dir_vep}/" && \
	cd "${cached_ref_dir_vep}/" && \
	tar xzf homo_sapiens_vep_101_GRCh38.tar.gz && \
	rm homo_sapiens_vep_101_GRCh38.tar.gz
	"""
}

// Depending on whether the reference files used for VEP annotation was pre-downloaded, set the input
// channel for the VEP annotation process
if( params.vep_ref_cached == "yes" ) {
	vep_ref_dir = vep_ref_dir_preDownloaded
}
else {
	vep_ref_dir = vep_ref_dir_fromProcess
}

// Combine all needed reference FASTA files into one channel for use in VEP annotation process
reference_genome_fasta_forAnnotation.combine( reference_genome_fasta_index_forAnnotation )
	.combine( reference_genome_fasta_dict_forAnnotation )
	.set{ reference_genome_bundle_forAnnotation }


/*

// VEP ~ annotate the final somatic VCFs using databases including Ensembl, GENCODE, RefSeq, PolyPhen, SIFT, dbSNP, COSMIC, etc.
process annotateSomaticVcf_vep {
	publishDir "${params.output_dir}/somatic/vepAnnotatedVcfs", mode: 'copy'
	tag "${vcf_id}"

	input:
	tuple path(high_quality_consensus_somatic_snv_vcf), path(cached_ref_dir_vep), path(reference_genome_fasta_forAnnotation), path(reference_genome_fasta_index_forAnnotation), path(reference_genome_fasta_dict_forAnnotation) from high_quality_consensus_snv_vcf_forAnnotation.combine(vep_ref_dir).combine(reference_genome_bundle_forAnnotation)

	output:
	path final_annotated_somatic_vcfs
	path annotation_summary

	when:
	params.varscan == "on" && params.mutect == "on" && params.caveman == "on" && params.strelka == "on"

	script:
	vcf_id = "${high_quality_consensus_somatic_snv_vcf}".replaceFirst(/\.vcf\.gz/, "")
	final_annotated_somatic_vcfs = "${vcf_id}.annotated.vcf.gz"
	annotation_summary = "${high_quality_consensus_somatic_snv_vcf}".replaceFirst(/\.vcf\.gz/, ".vep.summary.html")
	"""
	zgrep -E "^#|PASS" "${high_quality_consensus_somatic_snv_vcf}" > "${vcf_id}.passonly.vcf"

	vep \
	--offline \
	--cache \
	--dir "${cached_ref_dir_vep}" \
	--assembly GRCh38 \
	--fasta "${reference_genome_fasta_forAnnotation}" \
	--input_file "${vcf_id}.passonly.vcf" \
	--format vcf \
	--hgvs \
	--hgvsg \
	--protein \
	--symbol \
	--ccds \
	--canonical \
	--biotype \
	--sift b \
	--polyphen b \
	--stats_file "${annotation_summary}" \
	--output_file "${final_annotated_somatic_vcfs}" \
	--compress_output bgzip \
	--vcf
	"""
}

*/


// END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \\

