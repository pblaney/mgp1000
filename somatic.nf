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

		nextflow run somatic.nf -bg -resume --run_id batch1 --sample_sheet samplesheet.csv --email someperson@gmail.com --mutect_ref_vcf_concatenated no --annotsv_ref_cached no --vep_ref_cached no -profile somatic 

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
		                                      Default: input/preprocessedBams/
		--output_dir                   [str]  Directory that will hold all output files from the somatic variant analysis
		                                      Default: output/
		--email                        [str]  Email address to send workflow completion/stoppage notification
		--mutect_ref_vcf_concatenated  [str]  Indicates whether or not the gnomAD allele frequency reference VCF used for MuTect2 processes has
		                                      been concatenated, this will be done in a process of the pipeline if it has not, this does not
		                                      need to be done for every separate run after the first
		                                      Available: yes, no
		                                      Default: yes
		--battenberg_ref_cached        [str]  Indicates whether or not the reference files used for Battenberg have been downloaded/cached
		                                      locally, this will be done in a process of the pipeline if it has not, this does not need to be
		                                      done for every separate run after the first
		                                      Available: yes, no
		                                      Default: yes
		--annotsv_ref_cached           [str]  Indicates whether or not the AnnotSV reference files used for annotation have been downloaded/cached
		                                      locally, this will be done in a process of the pipeline if it has not, this does not need to be
		                                      done for every separate run after the first
		                                      Available: yes, no
		                                      Default: yes
		--vep_ref_cached               [str]  Indicates whether or not the VEP reference files used for annotation have been downloaded/cached
		                                      locally, this will be done in a process of the pipeline if it has not, this does not need to be
		                                      done for every separate run after the first
		                                      Available: yes, no
		                                      Default: yes
		--read_length                  [int]  Manually set the read length to be used for the Accucopy algorithm and mappability track for Control-FREEC
		                                      Available: 85, 100, 101, 150, 151 etc
		                                      Default: 100
		--cpus                         [int]  Globally set the number of cpus to be allocated for all processes that allow for multithreading
		                                      Available: 2, 4, 8, 16, etc.
		                                      Default: uniquly set for each process in nextflow.config to minimize resources needed
		--memory                       [str]  Globally set the amount of memory to be allocated for all processes, written as '##.GB' or '##.MB'
		                                      Available: 32.GB, 2400.MB, etc.
		                                      Default: uniquly set for each process in nextflow.config to minimize resources needed
		--queue_size                   [int]  Set max number of tasks the pipeline will handle in parallel
		                                      Available: 25, 50, 100, 150, etc.
		                                      Default: 100
		--executor                     [str]  Set the job executor for the run, this determines the system where the pipeline processes are run
		                                      and supervises their execution
		                                      Available: local, slurm
		                                      Default: slurm
		--help                        [flag]  Prints this message

	Toolbox Switches and Options:
		--telomerecat                  [str]  Indicates whether or not to use this tool
		                                      Available: off, on
		                                      Default: on
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
		--copycat                      [str]  Indicates whether or not to use this tool
		                                      Available: off, on
		                                      Default: on
		--battenberg                   [str]  Indicates whether or not to use this tool
		                                      Available: off, on
		                                      Default: on
		--battenberg_min_depth         [int]  Manually set the minimum read depth in the normal sample for SNP filtering in BAF calculations
		                                      Available: 3, 5, 7, 10
		                                      Default: 10
		--controlfreec                 [str]  Indicates whether or not to use this tool
		                                      Available: off, on
		                                      Default: on
		--controlfreec_bp_threshold  [float]  Manually set the breakpoint threshold value to be used for the Control-FREEC algorithm, this can be lowered
		                                      if the sample is expected to have large number of CNV segments or increased for the opposite assumption
		                                      Available: 0.6, 0.8, 1.2
		                                      Default: 0.8
		--controlfreec_ploidy          [int]  Manually set the ploidy value to be used for the Control-FREEC algorithm
		                                      Available: 2, 3, 4
		                                      Default: 2
		--sclust                       [str]  Indicates whether or not to use this tool
		                                      Available: off, on
		                                      Default: on
		--sclust_minp                [float]  Manually set the minimal expected ploidy to be used for the Sclust algorithm
		                                      Available: 1.5, 2.0
		                                      Default: 1.5
		--sclust_maxp                [float]  Manually set the maximal expected ploidy to be used for the Sclust algorithm
		                                      Available: 2.0, 3.5, 4.5, etc
		                                      Default: 4.5
		--sclust_mutclustering         [str]  Manually turn on or off the mutational clustering step of the Sclust process, this can be done if the process
		                                      cannot reach a solution for a given sample, this should only be used after attempts at lowering the lambda value
		                                      does not work, see --sclust_lambda parameter
		                                      Available: off, on
		                                      Default: on
		--sclust_lambda                [str]  Manually set the degree of smoothing for clustering mutations, increasing the value should resolve
		                                      issues with QP iterations related errors
		                                      Available: 1e-6, 1e-5
		                                      Default: 1e-7
		--accucopy                     [str]  Indicates whether or not to use this tool
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
		--igcaller                     [str]  Indicates whether or not to use this tool
		                                      Available: off, on
		                                      Default: on

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
params.email = null
params.mutect_ref_vcf_concatenated = "yes"
params.battenberg_ref_cached = "yes"
params.annotsv_ref_cached = "yes"
params.vep_ref_cached = "yes"
params.read_length = 100
params.telomerecat = "on"
params.telomerehunter = "on"
params.conpair = "on"
params.varscan = "on"
params.mutect = "on"
params.strelka = "on"
params.copycat = "on"
params.battenberg = "on"
params.controlfreec = "on"
params.sclust = "on"
params.manta = "on"
params.svaba = "on"
params.delly = "on"
params.igcaller = "on"
params.battenberg_min_depth = 10
params.controlfreec_bp_threshold = 0.8
params.controlfreec_ploidy = 2
params.sclust_minp = 1.5
params.sclust_maxp = 4.5
params.sclust_mutclustering = "on"
params.sclust_lambda = null
params.accucopy = "on"
params.cpus = null
params.memory = null
params.queue_size = 100
params.executor = 'slurm'
params.help = null

// Print help message if requested
if( params.help ) exit 0, helpMessage()

// Print preemptive error message if user-defined input/output directories does not exist
if( !file(params.input_dir).exists() ) exit 1, "The user-specified input directory does not exist in filesystem."

// Print preemptive error messages if required parameters are not set
if( params.run_id == null ) exit 1, "The run command issued does not have the '--run_id' parameter set. Please set the '--run_id' parameter to a unique identifier for the run."

if( params.sample_sheet == null ) exit 1, "The run command issued does not have the '--sample_sheet' parameter set. Please set the '--sample_sheet' parameter to the path of the normal/tumor pair sample sheet CSV."

// Print preemptive error message if Sclust is set while Mutect2 is not
if( params.sclust == "on" && params.mutect == "off" ) exit 1, "Sclust requires output from Mutect2 to run so both must be turned on"

// Print preemptive error message if Strelka is set while Manta is not
if( params.strelka == "on" && params.manta == "off" ) exit 1, "Strelka requires output from Manta to run so both must be turned on"

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
	       reference_genome_fasta_forControlFreecSamtoolsMpileup;
	       reference_genome_fasta_forControlFreecCalling;
	       reference_genome_fasta_forAccucopy;
	       reference_genome_fasta_forManta;
	       reference_genome_fasta_forStrelka;
	       reference_genome_fasta_forStrelkaBcftools;
	       reference_genome_fasta_forSvabaBcftools;
	       reference_genome_fasta_forDelly;
	       reference_genome_fasta_forIgCaller;
	       reference_genome_fasta_forConsensusSnvMpileup;
	       reference_genome_fasta_forConsensusIndelMpileup;
	       reference_genome_fasta_forConsensusSvFpFilter;
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
	       reference_genome_fasta_index_forControlFreecSamtoolsMpileup;
	       reference_genome_fasta_index_forControlFreecCalling;
	       reference_genome_fasta_index_forControlFreecConsensusPrep;
	       reference_genome_fasta_index_forSclustConsensusCnv;
	       reference_genome_fasta_index_forAccucopy;
	       reference_genome_fasta_index_forManta;
	       reference_genome_fasta_index_forStrelka;
	       reference_genome_fasta_index_forStrelkaBcftools;
	       reference_genome_fasta_index_forSvabaBcftools;
	       reference_genome_fasta_index_forDelly;
	       reference_genome_fasta_index_forIgCaller;
	       reference_genome_fasta_index_forConsensusSnvMpileup;
	       reference_genome_fasta_index_forConsensusIndelMpileup;
	       reference_genome_fasta_index_forConsensusSvFpFilter;
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
	       reference_genome_fasta_dict_forControlFreecSamtoolsMpileup;
	       reference_genome_fasta_dict_forControlFreecCalling;
	       reference_genome_fasta_dict_forAccucopy;
	       reference_genome_fasta_dict_forManta;
	       reference_genome_fasta_dict_forStrelka;
	       reference_genome_fasta_dict_forStrelkaBcftools;
	       reference_genome_fasta_dict_forSvaba;
	       reference_genome_fasta_dict_forSvabaBcftools;
	       reference_genome_fasta_dict_forDelly;
	       reference_genome_fasta_dict_forIgCaller;
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
	       gatk_bundle_wgs_bed_forStrelka }

Channel
	.fromPath( 'references/hg38/wgs_calling_regions_blacklist.0based.hg38.bed' )
	.into{ gatk_bundle_wgs_bed_blacklist_0based_forDelly;
	       gatk_bundle_wgs_bed_blacklist_0based_forSvaba }

Channel
	.fromList( ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6',
	            'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
	            'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18',
	            'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY',] )
	.into{ chromosome_list_forVarscanSamtoolsMpileup;
	       chromosome_list_forMutectCalling;
	       chromosome_list_forMutectPileup;
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

if( params.battenberg == "on" && params.battenberg_ref_cached == "yes" ) {
	Channel
		.fromPath( 'references/hg38/battenberg_reference/', checkIfExists: true )
		.ifEmpty{ error "The run command issued has the '--battenberg_ref_cached' parameter set to 'yes', however the directory does not exist. Please set the '--battenberg_ref_cached' parameter to 'no' and resubmit the run command. For more information, check the README or issue the command 'nextflow run somatic.nf --help'"}
		.set{ battenberg_ref_dir_preDownloaded }
}

Channel
	.fromPath( 'references/hg38/SnpGcCorrections.hg38.tsv' )
	.set{ snp_gc_corrections }

Channel
	.fromPath( 'references/hg38/Homo_sapiens_assembly38_autosome_sex_chroms', type: 'dir' )
	.set{ autosome_sex_chromosome_fasta_dir }

Channel
	.fromPath( 'references/hg38/Homo_sapiens_assembly38_autosome_sex_chrom_sizes.txt' )
	.into{ autosome_sex_chromosome_sizes_forControlFreec;
	       autosome_sex_chromosome_sizes_forCopycat }

Channel
	.fromPath( 'references/hg38/common_all_20180418.vcf.gz' )
	.set{ common_dbsnp_ref_vcf }

Channel
	.fromPath( 'references/hg38/common_all_20180418.vcf.gz.tbi' )
	.set{ common_dbsnp_ref_vcf_index }

Channel
	.fromPath( 'references/hg38/mappability_track_85m2.hg38.zip' )
	.set{ mappability_track_85kmer_zip }

Channel
	.fromPath( 'references/hg38/mappability_track_100m2.hg38.zip' )
	.set{ mappability_track_100kmer_zip }

Channel
	.fromPath( 'references/hg38/mappability_track_150m2.hg38.zip' )
	.set{ mappability_track_150kmer_zip }

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
	.fromPath( 'references/hg38/snp_sites_common_1000_genomes.hg38.gz' )
	.set{ common_1000G_snps_sites }

Channel
	.fromPath( 'references/hg38/snp_sites_common_1000_genomes.hg38.gz.tbi' )
	.set{ common_1000G_snps_sites_index }

Channel
	.fromPath( 'references/hg38/simple_and_centromeric_repeats.hg38.bed' )
	.into{ simple_and_centromeric_repeats_bed_forSvaba;
		   simple_and_centromeric_repeats_bed_forSnvBedFilter;
	       simple_and_centromeric_repeats_bed_forIndelBedFilter }

if( params.annotsv_ref_cached == "yes" ) {
     Channel
          .fromPath( 'references/hg38/annotations_human_annotsv_hg38/', type: 'dir', checkIfExists: true )
          .ifEmpty{ error "The run command issued has the '--annotsv_ref_cached' parameter set to 'yes', however the directory does not exist. Please set the '--annotsv_ref_cached' parameter to 'no' and resubmit the run command. For more information, check the README or issue the command 'nextflow run somatic.nf --help'"}
          .into{ annotsv_ref_dir_pre_downloaded_forSvAnnotation;
                 annotsv_ref_dir_pre_downloaded_forCnvAnnotation }
}

if( params.vep_ref_cached == "yes" ) {
	Channel
		.fromPath( 'references/hg38/homo_sapiens_vep_101_GRCh38/', type: 'dir', checkIfExists: true )
		.ifEmpty{ error "The run command issued has the '--vep_ref_cached' parameter set to 'yes', however the directory does not exist. Please set the '--vep_ref_cached' parameter to 'no' and resubmit the run command. For more information, check the README or issue the command 'nextflow run somatic.nf --help'"}
		.set{ vep_ref_dir_preDownloaded }
}

// #################################################### \\
// ~~~~~~~~~~~~~~~~ PIPELINE PROCESSES ~~~~~~~~~~~~~~~~ \\

log.info ''
log.info '######### Myeloma Genome Pipeline 1000 #########'
log.info '################################################'
log.info '~~~~~~~~~~~~~~~~~~~ SOMATIC ~~~~~~~~~~~~~~~~~~~'
log.info '################################################'
log.info ''
log.info "~~~ Launch Time ~~~		${workflowTimestamp}"
log.info ''
log.info "~~~ Input Directory ~~~		${params.input_dir}"
log.info ''
log.info "~~~ Output Directory ~~~	${params.output_dir}"
log.info ''
log.info "~~~ Run Report File ~~~		nextflow_report.${params.run_id}.html"
log.info ''
log.info "~~~ Read Length ~~~		${params.read_length}"
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
		   tumor_normal_pair_forTelomerecat;
		   tumor_normal_pair_forTelomereHunter;
		   tumor_normal_pair_forConpairPileup;
	       tumor_normal_pair_forVarscanSamtoolsMpileup;
	       tumor_normal_pair_forMutectCalling;
	       tumor_normal_pair_forMutectPileup;
	       tumor_normal_pair_forControlFreecSamtoolsMpileup;
	       tumor_normal_pair_forSclustBamprocess;
	       tumor_normal_pair_forAccucopy;
	       tumor_normal_pair_forManta;
	       tumor_normal_pair_forSvaba;
	       tumor_normal_pair_forDelly;
	       tumor_normal_pair_forIgCaller }

// Combine reference FASTA index and sex identification loci files into one channel for use in alleleCount process
reference_genome_fasta_index_forAlleleCount.combine( sex_identification_loci )
	.set{ ref_index_and_sex_ident_loci }

// alleleCount ~ determine the sex of each sample to use in downstream analyses
process identifySampleSex_allelecount {
	publishDir "${params.output_dir}/somatic/sexOfSamples", mode: 'copy', pattern: '*.{txt}'
	tag "${tumor_normal_sample_id}"

	input:
	tuple path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(reference_genome_fasta_index_forAlleleCount), path(sex_identification_loci) from tumor_normal_pair_forAlleleCount.combine(ref_index_and_sex_ident_loci)

	output:
	tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index) into bams_forVarscanBamReadcount
	tuple val(tumor_normal_sample_id), path(sample_sex) into sex_of_sample_forControlFreecCalling
	tuple val(tumor_normal_sample_id), path(tumor_bam) into bam_forCopycat
	tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(sample_sex) into bams_and_sex_of_sample_forBattenberg
	tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index) into bams_forConsensusSnvMpileup
	tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index) into bams_forConsensusIndelMpileup
	tuple val(tumor_normal_sample_id), path(sample_sex) into sex_of_sample_forConsensusCnvTransform
	tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index) into bam_forConsensusSvFpFilter
	tuple val(tumor_normal_sample_id), path(sample_sex) into allelecount_output_forConsensusMetadata

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

// Telomerecat bam2length ~  estimating the average telomere length
process telomereLengthEstimation_telomerecat {
    publishDir "${params.output_dir}/somatic/telomerecat", mode: 'copy', pattern: '*.{csv}'
    tag "${tumor_normal_sample_id}"

    input:
    tuple path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index) from tumor_normal_pair_forTelomerecat

    output:
    path normal_telomere_estimates
    path tumor_telomere_estimates

    when:
    params.telomerecat == "on"

    script:
    tumor_id = "${tumor_bam.baseName}".replaceFirst(/\..*$/, "")
    normal_id = "${normal_bam.baseName}".replaceFirst(/\..*$/, "")
    tumor_normal_sample_id = "${tumor_id}_vs_${normal_id}"
    normal_telomere_estimates = "${normal_id}.telomerecat.csv"
    tumor_telomere_estimates = "${tumor_id}.telomerecat.csv"
    """
    telomerecat bam2length \
    -p ${task.cpus} \
    -v 1 \
    --output "${tumor_telomere_estimates}" \
    "${tumor_bam}"

    telomerecat bam2length \
    -p ${task.cpus} \
    -v 1 \
    --output "${normal_telomere_estimates}" \
    "${normal_bam}"
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
	tag "${tumor_normal_sample_id}"

	input:
	tuple path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(reference_genome_fasta_forConpairPileup), path(reference_genome_fasta_index_forConpairPileup), path(reference_genome_fasta_dict_forConpairPileup) from tumor_normal_pair_forConpairPileup.combine(reference_genome_bundle_forConpairPileup)

	output:
	tuple val(tumor_normal_sample_id), path(tumor_pileup), path(normal_pileup) into bam_pileups_forConpair

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
	publishDir "${params.output_dir}/somatic/conpair", mode: 'copy', pattern: '*.{txt}'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(tumor_pileup), path(normal_pileup) from bam_pileups_forConpair

	output:
	tuple val(tumor_normal_sample_id), path(conpair_concordance_file), path(conpair_contamination_file) into conpair_output_forConsensusMetadata

	when:
	params.conpair == "on"
	
	script:
	conpair_concordance_file = "${tumor_normal_sample_id}.conpair.concordance.txt"
	conpair_contamination_file = "${tumor_normal_sample_id}.conpair.contamination.txt"
	hg38_ref_genome_markers = "/data/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover"
	"""
	\${CONPAIR_DIR}/scripts/verify_concordance.py \
	--min_cov 10 \
	--min_mapping_quality 10 \
	--min_base_quality 20 \
	--tumor_pileup "${tumor_pileup}" \
	--normal_pileup "${normal_pileup}" \
	--outfile "${conpair_concordance_file}" \
	--markers \${CONPAIR_DIR}"${hg38_ref_genome_markers}.txt"

	\${CONPAIR_DIR}/scripts/estimate_tumor_normal_contamination.py \
	--grid 0.01 \
	--min_mapping_quality 10 \
	--tumor_pileup "${tumor_pileup}" \
	--normal_pileup "${normal_pileup}" \
	--outfile "${conpair_contamination_file}" \
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
	tag "${tumor_normal_sample_id} C=${chromosome}"

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


// ~~~~~~~~~~~~~~~~ Mutect2 ~~~~~~~~~~~~~~~~ \\
// START

// BCFtools Concat ~ prepare the gnomAD allele frequency reference VCF for Mutect2 process, if needed
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
// channel for the for Mutect2 process
if( params.mutect_ref_vcf_concatenated == "yes" && params.mutect == "on") {
	mutect_gnomad_ref_vcf = mutect_gnomad_ref_vcf_preBuilt.combine( mutect_gnomad_ref_vcf_index_preBuilt )
}
else {
	mutect_gnomad_ref_vcf = mutect_gnomad_ref_vcf_fromProcess
}

// Combine all reference FASTA files, WGS BED file, and resource VCFs into one channel for use in Mutect2 calling process
reference_genome_fasta_forMutectCalling.combine( reference_genome_fasta_index_forMutectCalling )
	.combine( reference_genome_fasta_dict_forMutectCalling )
	.combine( gatk_bundle_wgs_bed_forMutectCalling )
	.set{ reference_genome_bundle_and_bed_forMutectCalling }

reference_genome_bundle_and_bed_forMutectCalling.combine( mutect_gnomad_ref_vcf )
	.combine( panel_of_normals_1000G )
	.combine( panel_of_normals_1000G_index )
	.set{ reference_genome_bed_and_vcfs_forMutectCalling } 

// GATK Mutect2 ~ call somatic SNVs and indels via local assembly of haplotypes
process snvAndIndelCalling_gatk {
	tag "${tumor_normal_sample_id} C=${chromosome} "

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

// GATK SortVcfs ~ merge all per chromosome Mutect2 VCFs
process mergeAndSortMutect2Vcfs_gatk {
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

// GATK MergeMutectStats ~ merge the per chromosome stats output of Mutect2 variant calling
process mergeMutect2StatsForFiltering_gatk {
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

// Combine WGS BED file and resource VCFs into one channel for use in Mutect2 pileup process
gatk_bundle_wgs_bed_forMutectPileup.combine( exac_common_sites_ref_vcf )
	.combine( exac_common_sites_ref_vcf_index )
	.set{ bed_and_resources_vcfs_forMutectPileup }

// GATK GetPileupSummaries ~ tabulates pileup metrics for inferring contamination
process pileupSummariesForMutect2Contamination_gatk {
	tag "${tumor_normal_sample_id} C=${chromosome}"

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
	publishDir "${params.output_dir}/somatic/mutect", mode: 'copy', pattern: '*.{txt}'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(tumor_pileup), path(normal_pileup) from tumor_pileups_forMutectContamination.join(normal_pileups_forMutectContamination)

	output:
	tuple val(tumor_normal_sample_id), path(mutect_contamination_file) into contamination_file_forMutectFilter, mutect_output_forConsensusMetadata

	when:
	params.mutect == "on"

	script:
	mutect_contamination_file = "${tumor_normal_sample_id}.mutect.contamination.txt" 
	"""
	gatk CalculateContamination \
	--java-options "-Xmx${task.memory.toGiga()}G -Djava.io.tmpdir=." \
	--verbosity ERROR \
	--tmp-dir . \
	--input "${tumor_pileup}" \
	--matched-normal "${normal_pileup}" \
	--output "${mutect_contamination_file}"
	"""
}

// Combine all reference FASTA files and input VCF, stats file, and contamination table into one channel for use in Mutect2 filtering process
reference_genome_fasta_forMutectFilter.combine( reference_genome_fasta_index_forMutectFilter )
	.combine( reference_genome_fasta_dict_forMutectFilter )
	.set{ reference_genome_bundle_forMutectFilter }

merged_raw_vcfs_forMutectFilter.join( merged_mutect_stats_file_forMutectFilter )
	.join( contamination_file_forMutectFilter )
	.set{ input_vcf_stats_and_contamination_forMutectFilter }

// GATK FilterMutectCalls ~ filter somatic SNVs and indels called by Mutect2
process mutect2VariantFiltration_gatk {
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(merged_raw_vcf), path(merged_raw_vcf_index), path(merged_mutect_stats_file), path(mutect_contamination_file), path(reference_genome_fasta_forMutectFilter), path(reference_genome_fasta_index_forMutectFilter), path(reference_genome_fasta_dict_forMutectFilter) from input_vcf_stats_and_contamination_forMutectFilter.combine(reference_genome_bundle_forMutectFilter)

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
	--contamination-table "${mutect_contamination_file}" \
	--output "${filtered_vcf}" \
	--filtering-stats "${filter_stats_file}"
	"""
}

// Combine all needed reference FASTA files into one channel for use in Mutect2 / BCFtools Norm process
reference_genome_fasta_forMutectBcftools.combine( reference_genome_fasta_index_forMutectBcftools )
	.combine( reference_genome_fasta_dict_forMutectBcftools )
	.set{ reference_genome_bundle_forMutectBcftools }

// BCFtools Norm ~ split multiallelic sites into multiple rows then left-align and normalize indels
process splitMultiallelicAndLeftNormalizeMutect2Vcf_bcftools {
	publishDir "${params.output_dir}/somatic/mutect", mode: 'copy', pattern: '*.{txt}'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(filtered_vcf), path(filtered_vcf_index), path(reference_genome_fasta_forMutectBcftools), path(reference_genome_fasta_index_forMutectBcftools), path(reference_genome_fasta_dict_forMutectBcftools) from filtered_vcf_forMutectBcftools.combine(reference_genome_bundle_forMutectBcftools)

	output:
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


// ~~~~~~~~~~~~~~~~ Copycat ~~~~~~~~~~~~~~~~ \\
// START

// Copycat ~ capture and bin the read coverage across a genome for CNV and SV support
process binReadCoverage_copycat {
    publishDir "${params.output_dir}/somatic/copycat", mode: 'copy', pattern: '*.{csv.gz,seg}'
    tag "${tumor_normal_sample_id}"

    input:
    tuple val(tumor_normal_sample_id), path(tumor_bam), path(autosome_sex_chromosome_sizes) from bam_forCopycat.combine(autosome_sex_chromosome_sizes_forCopycat)
     
    output:
    path tumor_copycat_coverage_csv
    path tumor_copycat_coverage_seg

    when:
    params.copycat == "on"

    script:
    tumor_id = "${tumor_bam.baseName}".replaceFirst(/\..*$/, "")
    tumor_copycat_coverage_csv = "${tumor_id}.coverage.10kb.csv.gz"
    tumor_copycat_coverage_seg = "${tumor_id}.coverage.10kb.igv.seg"
    """
    copycat.sh \
    "${tumor_bam}" \
    "${autosome_sex_chromosome_sizes}" \
    "${tumor_id}"

    mv "${tumor_id}.coverage.10kb.for_IGV.seg" "${tumor_copycat_coverage_seg}"
    """
}

// END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \\


// ~~~~~~~~~~~~~~~ Battenberg ~~~~~~~~~~~~~~ \\
// START

// Battenberg ~ download the reference files needed to run Battenberg
process downloadBattenbergReferences_battenberg {
  	publishDir "references/hg38", mode: 'copy'

  	output:
  	path battenberg_references into battenberg_ref_dir_fromProcess

  	when:
  	params.battenberg == "on" && params.battenberg_ref_cached == "no"

  	script:
  	battenberg_references = "battenberg_reference"
  	"""
  	mkdir -p "${battenberg_references}"
  	cd "${battenberg_references}/"
  	mkdir -p GC_correction_hg38/
  	mkdir -p RT_correction_hg38/

  	battenberg_reference_downloader.sh
  	"""
}

// Depending on whether the reference files used for Battenberg was pre-downloaded, set the input
// channel for the Battenberg process
if( params.battenberg_ref_cached == "yes" ) {
	battenberg_ref_dir = battenberg_ref_dir_preDownloaded
}
else {
	battenberg_ref_dir = battenberg_ref_dir_fromProcess
}

// Battenberg ~ whole genome sequencing subclonal copy number caller
process cnvCalling_battenberg {
  	publishDir "${params.output_dir}/somatic/battenberg", mode: 'copy'
  	tag "${tumor_normal_sample_id}"

  	input:
  	tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(sample_sex), path(battenberg_references) from bams_and_sex_of_sample_forBattenberg.combine(battenberg_ref_dir)

  	output:
  	tuple val(tumor_normal_sample_id), path(battenberg_cnv_profile) into final_battenberg_cnv_profile_forConsensusPrep
  	tuple path(battenberg_rho_and_psi), path(battenberg_purity_and_ploidy), path(battenberg_cnv_profile_png)
  	path "${output_dir}/*.tumour.png"
  	path "${output_dir}/*.germline.png"
  	path "${output_dir}/*_distance.png"
  	path "${output_dir}/*_coverage.png"
  	path "${output_dir}/*_alleleratio.png"
  	path "${output_dir}/*_BattenbergProfile_subclones.png"
  	path "${output_dir}/*chr*_heterozygousData.png"
  	path "${output_dir}/*_nonroundedprofile.png"
  	path "${output_dir}/*_segment_chr*.png"
  	path "${output_dir}/*_GCwindowCorrelations_*Correction.txt"
  	path "${output_dir}/*_refit_suggestion.txt"
  	path "${output_dir}/*_segment_masking_details.txt"
  	path "${output_dir}/*_subclones_alternatives.txt"
  	path "${output_dir}/*_subclones_chr*.png"
  	path "${output_dir}/*_totalcn_chrom_plot.png"
  	path "${output_dir}/*_copynumberprofile.png"

  	when:
  	params.battenberg == "on"

  	script:
  	tumor_id = "${tumor_bam.baseName}".replaceFirst(/\..*$/, "")
  	normal_id = "${normal_bam.baseName}".replaceFirst(/\..*$/, "")
  	output_dir = "${tumor_normal_sample_id}_results"
  	battenberg_cnv_profile = "${tumor_normal_sample_id}.battenberg.cnv.txt"
  	battenberg_rho_and_psi = "${tumor_normal_sample_id}.battenberg.rho.psi.txt"
  	battenberg_purity_and_ploidy = "${tumor_normal_sample_id}.battenberg.purity.ploidy.txt"
  	battenberg_cnv_profile_png = "${tumor_normal_sample_id}.battenberg.cnv.png"
  	"""
  	sex=\$(cut -d ' ' -f 1 "${sample_sex}")

  	battenberg_executor.sh \
  	"${tumor_id}" \
  	"${normal_id}" \
  	"${tumor_bam}" \
  	"${normal_bam}" \
  	"\${sex}" \
  	"${output_dir}" \
  	"${task.cpus}" \
  	"${params.battenberg_min_depth}"

  	cp "${output_dir}/${tumor_id}_subclones.txt" "${battenberg_cnv_profile}"
  	cp "${output_dir}/${tumor_id}_rho_and_psi.txt" "${battenberg_rho_and_psi}"
  	cp "${output_dir}/${tumor_id}_purity_ploidy.txt" "${battenberg_purity_and_ploidy}"
  	cp "${output_dir}/${tumor_id}_BattenbergProfile_average.png" "${battenberg_cnv_profile_png}"
  	"""
}

// Battenberg Consensus CNV Prep ~ extract and prepare CNV output for consensus
process consensusCnvPrep_battenberg {
    tag "${tumor_normal_sample_id}"

    input:
    tuple val(tumor_normal_sample_id), path(battenberg_cnv_profile) from final_battenberg_cnv_profile_forConsensusPrep

    output:
    tuple val(tumor_normal_sample_id), path(battenberg_somatic_cnv_bed), path(battenberg_somatic_alleles_bed) into final_battenberg_cnv_profile_forConsensus
    tuple val(tumor_normal_sample_id), path(battenberg_subclones_file) into battenberg_subclones_forConsensusSubclones

    when:
    params.battenberg == "on"

    script:
    battenberg_somatic_cnv_bed = "${tumor_normal_sample_id}.battenberg.somatic.cnv.bed"
    battenberg_somatic_alleles_bed = "${tumor_normal_sample_id}.battenberg.somatic.alleles.bed"
    battenberg_subclones_file = "${tumor_normal_sample_id}.battenberg.subclones.txt"
    """
    battenberg_clonal_and_subclonal_extractor.py \
    <(grep -v 'startpos' "${battenberg_cnv_profile}" | cut -f 1-3,8-13) \
    "${battenberg_somatic_cnv_bed}" \
    "${battenberg_somatic_alleles_bed}" \
    "${battenberg_subclones_file}"
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
	tag "${tumor_normal_sample_id} C=${chromosome}"

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

// Set the mappability track for Control-FREEC based on user-defined read length
if( params.read_length < 100 ) {
  	mappability_track_zip = mappability_track_85kmer_zip
} else if( params.read_length < 150 ) {
  	mappability_track_zip = mappability_track_100kmer_zip
} else {
  	mappability_track_zip = mappability_track_150kmer_zip
}

// Combine mpileup input files with the sample sex identity then all reference files into one channel for use in Control-FREEC
tumor_normal_pileups_forControlFreecCalling.join( sex_of_sample_forControlFreecCalling )
	.set{ tumor_normal_pileups_and_sex_ident }

reference_genome_fasta_forControlFreecCalling.combine( reference_genome_fasta_index_forControlFreecCalling )
	.combine( reference_genome_fasta_dict_forControlFreecCalling )
	.set{ reference_genome_bundle_forControlFreecCalling }

reference_genome_bundle_forControlFreecCalling.combine( autosome_sex_chromosome_fasta_dir )
	.combine( autosome_sex_chromosome_sizes_forControlFreec )
	.combine( common_dbsnp_ref_vcf )
	.combine( common_dbsnp_ref_vcf_index )
	.combine( mappability_track_zip )
	.set{ reference_data_bundle_forControlFreec }

// Control-FREEC ~ detection of copy-number changes and allelic imbalances
process cnvCalling_controlfreec {
	publishDir "${params.output_dir}/somatic/controlFreec", mode: 'copy', pattern: '*.{txt,cnv,bedGraph}'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(tumor_pileup), path(normal_pileup), path(sex_of_sample_forControlFreecCalling), path(reference_genome_fasta_forControlFreecCalling), path(reference_genome_fasta_index_forControlFreecCalling), path(reference_genome_fasta_dict_forControlFreecCalling), path(autosome_sex_chromosome_fasta_dir), path(autosome_sex_chromosome_sizes), path(common_dbsnp_ref_vcf), path(common_dbsnp_ref_vcf_index), path(mappability_track_zip) from tumor_normal_pileups_and_sex_ident.combine(reference_data_bundle_forControlFreec)

	output:
	tuple val(tumor_normal_sample_id), path(cnv_profile_raw), path(cnv_ratio_file), path(baf_file), path(control_freec_bedgraph) into cnv_calling_files_forControlFreecPostProcessing
	path control_freec_config_file
	tuple val(tumor_normal_sample_id), path(control_freec_subclones_file) into control_freec_subclones_forConsensusSubclones
	tuple val(tumor_normal_sample_id), path(control_freec_run_info) into control_freec_output_forConsensusMetadata

	when:
	params.controlfreec == "on"

	script:
	breakpoint_threshold = "${params.controlfreec_bp_threshold}"
	control_freec_ploidy = "${params.controlfreec_ploidy}"
	control_freec_config_file = "${tumor_normal_sample_id}.controlfreec.config.txt"
	control_freec_run_info = "${tumor_normal_sample_id}.controlfreec.runinfo.txt"
	cnv_profile_raw = "${tumor_normal_sample_id}.controlfreec.raw.cnv"
	cnv_ratio_file = "${tumor_normal_sample_id}.controlfreec.ratio.txt"
	control_freec_subclones_file = "${tumor_normal_sample_id}.controlfreec.subclones.txt"
	baf_file = "${tumor_normal_sample_id}.controlfreec.baf.txt"
	control_freec_bedgraph = "${tumor_normal_sample_id}.controlfreec.ratio.bedGraph"
	"""
	unzip -q "${mappability_track_zip}"
	gemFile=\$(ls -1 *.gem)
	sex=\$(cut -d ' ' -f 2 "${sex_of_sample_forControlFreecCalling}")

	touch "${control_freec_config_file}"
	echo "[general]" >> "${control_freec_config_file}"
	echo "BedGraphOutput = TRUE" >> "${control_freec_config_file}"
	echo "breakPointThreshold = ${breakpoint_threshold}" >> "${control_freec_config_file}"
	echo "breakPointType = 2" >> "${control_freec_config_file}"
	echo "chrFiles = \${PWD}/${autosome_sex_chromosome_fasta_dir}" >> "${control_freec_config_file}"
	echo "chrLenFile = \${PWD}/${autosome_sex_chromosome_sizes}" >> "${control_freec_config_file}"
	echo "contaminationAdjustment = TRUE" >> "${control_freec_config_file}"
	echo "forceGCcontentNormalization = 1" >> "${control_freec_config_file}"
	echo "gemMappabilityFile = \${PWD}/\${gemFile}" >> "${control_freec_config_file}"
	echo "minimalSubclonePresence = 20" >> "${control_freec_config_file}"
	echo "maxThreads = ${task.cpus}" >> "${control_freec_config_file}"
	echo "ploidy = ${control_freec_ploidy}" >> "${control_freec_config_file}"
	echo "sex = \${sex}" >> "${control_freec_config_file}"
	echo "uniqueMatch = TRUE" >> "${control_freec_config_file}"
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

	mv "${tumor_pileup}_info.txt" "${control_freec_run_info}"
	mv "${tumor_pileup}_CNVs" "${cnv_profile_raw}"
	mv "${tumor_pileup}_ratio.txt" "${cnv_ratio_file}"
	mv "${tumor_pileup}_subclones.txt" "${control_freec_subclones_file}"
	mv "${tumor_pileup}_BAF.txt" "${baf_file}"
	mv "${tumor_pileup}_ratio.BedGraph" "${control_freec_bedgraph}"
	"""
}

// Control-FREEC ~ post-processing of CNV predictions for significance, visualization, and format compatibility
process cnvPredictionPostProcessing_controlfreec {
	publishDir "${params.output_dir}/somatic/controlFreec", mode: 'copy', pattern: '*.{txt,bed,png}'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(cnv_profile_raw), path(cnv_ratio_file), path(baf_file), path(control_freec_bedgraph) from cnv_calling_files_forControlFreecPostProcessing

	output:
	tuple val(tumor_normal_sample_id), path(control_freec_bedgraph), path(control_freec_cnv_profile_final) into final_control_freec_cnv_profile_forConsensusPrep
	path ratio_graph_png
	path ratio_log2_graph_png
	path baf_graph_png

	when:
	params.controlfreec == "on"

	script:
	control_freec_cnv_profile_final = "${tumor_normal_sample_id}.controlfreec.cnv.txt"
	control_freec_cnv_ratio_bed_file = "${tumor_normal_sample_id}.controlfreec.ratio.bed"
	ratio_graph_png = "${tumor_normal_sample_id}.controlfreec.ratio.png"
	ratio_log2_graph_png = "${tumor_normal_sample_id}.controlfreec.ratio.log2.png"
	baf_graph_png = "${tumor_normal_sample_id}.controlfreec.baf.png"
	"""
	cat \${CONTROLFREEC_DIR}/scripts/assess_significance.R | R --slave --args "${cnv_profile_raw}" "${cnv_ratio_file}"
	mv "${cnv_profile_raw}.p.value.txt" "${control_freec_cnv_profile_final}"

	cat \${CONTROLFREEC_DIR}/scripts/makeGraph.R | R --slave --args 2 "${cnv_ratio_file}" "${baf_file}"
	mv "${cnv_ratio_file}.png" "${ratio_graph_png}"
	mv "${cnv_ratio_file}.log2.png" "${ratio_log2_graph_png}"
	mv "${baf_file}.png" "${baf_graph_png}"

	perl \${CONTROLFREEC_DIR}/scripts/freec2bed.pl -f "${cnv_ratio_file}" > "${control_freec_cnv_ratio_bed_file}"
	"""
}

// Control-FREEC Consensus CNV Prep ~ extract and prepare CNV output for consensus
process consensusCnvPrep_controlfreec {
  	tag "${tumor_normal_sample_id}"

  	input:
  	tuple val(tumor_normal_sample_id), path(control_freec_bedgraph), path(control_freec_cnv_profile_final), path(reference_genome_fasta_index_forControlFreecConsensusPrep) from final_control_freec_cnv_profile_forConsensusPrep.combine(reference_genome_fasta_index_forControlFreecConsensusPrep)

  	output:
  	tuple val(tumor_normal_sample_id), path(control_freec_somatic_cnv_bed), path(control_freec_somatic_alleles_bed) into final_controlfreec_cnv_profile_forConsensus

  	when:
  	params.controlfreec == "on"

  	script:
  	control_freec_somatic_cnv_bed = "${tumor_normal_sample_id}.controlfreec.somatic.cnv.bed"
  	control_freec_somatic_alleles_bed = "${tumor_normal_sample_id}.controlfreec.somatic.alleles.bed"
  	"""
  	control_freec_cnv_and_allele_preparer.sh \
  	"${tumor_normal_sample_id}" \
  	"${control_freec_cnv_profile_final}" \
  	"${control_freec_bedgraph}" \
  	"${reference_genome_fasta_index_forControlFreecConsensusPrep}"

  	control_freec_segment_refiner.py \
  	"${tumor_normal_sample_id}.controlfreec.complete.cnv.alleles.merged.bed" \
  	"${control_freec_somatic_alleles_bed}" \
  	"${control_freec_somatic_cnv_bed}"
  	"""
}

// END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \\


// ~~~~~~~~~~~~~~~~~ Sclust ~~~~~~~~~~~~~~~~~ \\
// START

// Sclust bamprocess ~ extract the read ratio and SNP information per chromosome
process bamprocessPerChromosome_sclust {
	tag "${tumor_normal_sample_id} C=${chromosome}"

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
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(final_mutect_vcf), path(final_mutect_vcf_index) from mutect_vcf_forSclust

	output:
	tuple val(tumor_normal_sample_id), path(mutations_vcf) into vcf_forSclustCn, vcf_forSclustClustering

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

// Sclust cn ~ perform copy-number analysis and estimation of tumor purity
process cnvCalling_sclust {
	publishDir "${params.output_dir}/somatic/sclust", mode: 'copy', pattern: '*.{txt,pdf}'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(read_count_file), path(common_snp_count_file), path(mutations_vcf) from read_count_and_snp_count_files.join(vcf_forSclustCn)

	output:
	tuple val(tumor_normal_sample_id), path(read_count_file), path(common_snp_count_file), path(sclust_allelic_states_file), path(sclust_subclones_file), path(sclust_cnv_summary_file), path(mutations_exp_af_file), path(sclust_cnv_segments_file) into sclust_cn_output_forClustering
	tuple val(tumor_normal_sample_id), path(sclust_allelic_states_file) into final_sclust_cnv_profile_forConsensusPrep
	tuple val(tumor_normal_sample_id), path(sclust_subclones_file) into sclust_subclones_forConsensusSubclones
	tuple val(tumor_normal_sample_id), path(sclust_cnv_summary_file) into sclust_output_forConsensusMetadata
	path mutations_exp_af_file
	path sclust_cnv_profile_pdf
	path sclust_cnv_segments_file

	when:
	params.sclust == "on" && params.mutect == "on"

	script:
	sclust_min_ploidy = "${params.sclust_minp}"
	sclust_max_ploidy = "${params.sclust_maxp}"
	sclust_lambda = params.sclust_lambda ? "-lambda ${params.sclust_lambda}" : ""
	sclust_allelic_states_file = "${tumor_normal_sample_id}.sclust.allelicstates.txt"
	sclust_cnv_profile_pdf = "${tumor_normal_sample_id}.sclust.profile.pdf"
	sclust_cnv_summary_file = "${tumor_normal_sample_id}.sclust.cnvsummary.txt"
	sclust_cnv_segments_file = "${tumor_normal_sample_id}.sclust.cnvsegments.txt"
	mutations_exp_af_file = "${tumor_normal_sample_id}_muts_expAF.txt"
	sclust_subclones_file = "${tumor_normal_sample_id}.sclust.subclones.txt"
	"""
	gunzip -f "${mutations_vcf}"

	Sclust cn \
	-ns 1000 \
	-rc "${read_count_file}" \
	-snp "${common_snp_count_file}" \
	-vcf "${tumor_normal_sample_id}.sclust.final.vcf" \
	-o "${tumor_normal_sample_id}" \
	-minp "${sclust_min_ploidy}" \
	-maxp "${sclust_max_ploidy}"

	mv "${tumor_normal_sample_id}_allelic_states.txt" "${sclust_allelic_states_file}"
	mv "${tumor_normal_sample_id}_cn_profile.pdf" "${sclust_cnv_profile_pdf}"
	mv "${tumor_normal_sample_id}_cn_summary.txt" "${sclust_cnv_summary_file}"
	mv "${tumor_normal_sample_id}_iCN.seg" "${sclust_cnv_segments_file}"
	mv "${tumor_normal_sample_id}_subclonal_cn.txt" "${sclust_subclones_file}"
	"""
}

// Sclust cluster ~ perform mutational clustering based on copy-number analysis and estimation of tumor purity
process mutationalClustering_sclust {
	publishDir "${params.output_dir}/somatic/sclust", mode: 'copy', pattern: '*.{txt,pdf}'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(read_count_file), path(common_snp_count_file), path(sclust_allelic_states_file), path(sclust_subclones_file), path(sclust_cnv_summary_file), path(mutations_exp_af_file), path(sclust_cnv_segments_file), path(mutations_vcf) from sclust_cn_output_forClustering.join(vcf_forSclustClustering)

	output:
	path mutation_clusters_file
	path mutation_clusters_pdf
	path cluster_assignment_file

	when:
	params.sclust == "on" && params.mutect == "on" && params.sclust_mutclustering == "on"

	script:
	sclust_lambda = params.sclust_lambda ? "-lambda ${params.sclust_lambda}" : ""
	mutation_clusters_file = "${tumor_normal_sample_id}.sclust.mutclusters.txt"
	mutation_clusters_pdf = "${tumor_normal_sample_id}.sclust.mutclusters.pdf"
	cluster_assignment_file = "${tumor_normal_sample_id}.sclust.clusterassignments.txt"
	"""
	gunzip -f "${mutations_vcf}"
	
	Sclust cluster \
	-i "${tumor_normal_sample_id}" \
	${sclust_lambda}

	mv "${tumor_normal_sample_id}_mclusters.txt" "${mutation_clusters_file}"
	mv "${tumor_normal_sample_id}_mcluster.pdf" "${mutation_clusters_pdf}"
	mv "${tumor_normal_sample_id}_cluster_assignments.txt" "${cluster_assignment_file}"
	"""
}

// Sclust Consensus CNV Prep ~ extract and prepare CNV output for consensus
process consensusCnvPrep_sclust {
  	tag "${tumor_normal_sample_id}"

  	input:
  	tuple val(tumor_normal_sample_id), path(sclust_allelic_states_file), path(reference_genome_fasta_index_forSclustConsensusCnv) from final_sclust_cnv_profile_forConsensusPrep.combine(reference_genome_fasta_index_forSclustConsensusCnv)

  	output:
  	tuple val(tumor_normal_sample_id), path(sclust_somatic_cnv_bed), path(sclust_somatic_alleles_bed) into final_sclust_cnv_profile_forConsensus

  	when:
  	params.sclust == "on" && params.mutect == "on"

  	script:
  	sclust_somatic_cnv_bed = "${tumor_normal_sample_id}.sclust.somatic.cnv.bed"
  	sclust_somatic_alleles_bed = "${tumor_normal_sample_id}.sclust.somatic.alleles.bed"
  	"""
  	# total copy number per segment
  	sclust_profile_gap_filler.py \
  	<(grep -v 'Sample' "${sclust_allelic_states_file}" | cut -f 2-4,6) \
  	<(head -n 24 "${reference_genome_fasta_index_forSclustConsensusCnv}" | cut -f 2) \
  	| \
  	sort -k1,1V -k2,2n > "${sclust_somatic_cnv_bed}"

  	# major/minor alleles per segment
  	sclust_profile_gap_filler.py \
  	<(grep -v 'Sample' "${sclust_allelic_states_file}" | awk 'BEGIN {OFS="\t"} {print \$2,\$3,\$4,\$7"/"\$8}') \
  	<(head -n 24 "${reference_genome_fasta_index_forSclustConsensusCnv}" | cut -f 2) \
  	| \
  	sort -k1,1V -k2,2n > "${sclust_somatic_alleles_bed}"
  	"""
}

// END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \\


// ~~~~~~~~~~~~~~~~ Accucopy ~~~~~~~~~~~~~~~~ \\
// START

// Combine reference FASTA and 1000 Genomes common SNPs file into one channel for use in Accucopy process
reference_genome_fasta_forAccucopy.combine( reference_genome_fasta_index_forAccucopy )
	.combine( reference_genome_fasta_dict_forAccucopy )
	.combine( common_1000G_snps_sites )
	.combine( common_1000G_snps_sites_index )
	.set{ ref_genome_and_snp_sites_forAccucopy }

// Accucopy ~ inference of allele-specific copy number alterations
process cnvCalling_accucopy {
	publishDir "${params.output_dir}/somatic/accucopy", mode: 'copy'
  	tag "${tumor_normal_sample_id}"

  	input:
  	tuple path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(reference_genome_fasta_forAccucopy), path(reference_genome_fasta_index_forAccucopy), path(reference_genome_fasta_dict_forAccucopy), path(common_1000G_snps_sites), path(common_1000G_snps_sites_index) from tumor_normal_pair_forAccucopy.combine(ref_genome_and_snp_sites_forAccucopy)

  	output:
  	tuple val(tumor_normal_sample_id), path(accucopy_cnv_profile) into accucopy_cnv_profile_forConsensusPrep
  	path "${tumor_normal_sample_id}/${accucopy_config_file}"
  	path accucopy_run_summary
  	path accucopy_detailed_run_status_log
  	path accucopy_cnv_png
  	path accucopy_tre_png
  	path accucopy_het_snps
  	path accucopy_segments
  	path "${tumor_normal_sample_id}/chr*.ratio.*.csv.gz"

  	when:
  	params.accucopy == "on"

  	script:
  	tumor_id = "${tumor_bam.baseName}".replaceFirst(/\..*$/, "")
	normal_id = "${normal_bam.baseName}".replaceFirst(/\..*$/, "")
	tumor_normal_sample_id = "${tumor_id}_vs_${normal_id}"
  	accucopy_config_file = "${tumor_normal_sample_id}.accucopy.config"
  	accucopy_run_summary = "${tumor_normal_sample_id}/${tumor_normal_sample_id}.accucopy.summary.txt"
	accucopy_detailed_run_status_log = "${tumor_normal_sample_id}/${tumor_normal_sample_id}.accucopy.status.txt"
	accucopy_cnv_profile = "${tumor_normal_sample_id}/${tumor_normal_sample_id}.accucopy.cnv.txt"
	accucopy_cnv_png = "${tumor_normal_sample_id}/${tumor_normal_sample_id}.accucopy.cnv.png"
	accucopy_tre_png = "${tumor_normal_sample_id}/${tumor_normal_sample_id}.accucopy.tre.png"
	accucopy_het_snps = "${tumor_normal_sample_id}/${tumor_normal_sample_id}.accucopy.het.snps.txt.gz"
	accucopy_segments = "${tumor_normal_sample_id}/${tumor_normal_sample_id}.accucopy.segments.txt.gz"
  	"""
  	if [[ ! "${tumor_bam_index}" =~ .bam.bai\$ ]]; then
    	cp "${tumor_bam_index}" "${tumor_bam}.bai"
  	fi
  	if [[ ! "${normal_bam_index}" =~ .bam.bai\$ ]]; then
    	cp "${normal_bam_index}" "${normal_bam}.bai"
  	fi

  	mkdir -p ref_dir/
  	cp "${reference_genome_fasta_forAccucopy}" ref_dir/genome.fa
  	cp "${reference_genome_fasta_index_forAccucopy}" ref_dir/genome.fa.fai
  	cp "${reference_genome_fasta_dict_forAccucopy}" ref_dir/genome.dict
  	cp "${common_1000G_snps_sites}" ref_dir/snp_sites.gz
  	cp "${common_1000G_snps_sites_index}" ref_dir/snp_sites.gz.tbi

  	touch "${accucopy_config_file}"
  	echo "read_length	${params.read_length}" >> "${accucopy_config_file}"
  	echo "window_size	500" >> "${accucopy_config_file}"
  	echo "reference_folder_path	\${PWD}/ref_dir" >> "${accucopy_config_file}"
  	echo "samtools_path	/usr/local/bin/samtools" >> "${accucopy_config_file}"
  	echo "caller_path	/usr/local/strelka" >> "${accucopy_config_file}"
  	echo "accucopy_path	/usr/local/Accucopy" >> "${accucopy_config_file}"

  	\${ACCUCOPY_DIR}/main.py \
  	--configure_filepath ${accucopy_config_file} \
  	--tumor_bam "${tumor_bam}" \
  	--normal_bam "${normal_bam}" \
  	--output_dir "${tumor_normal_sample_id}/" \
  	--nCores "${task.cpus}" \
  	--debug 1

  	cp "${accucopy_config_file}" "${tumor_normal_sample_id}/${accucopy_config_file}"
  	mv "${tumor_normal_sample_id}/infer.out.tsv" "${accucopy_run_summary}"
	mv "${tumor_normal_sample_id}/infer.status.txt" "${accucopy_detailed_run_status_log}"
	mv "${tumor_normal_sample_id}/cnv.output.tsv" "${accucopy_cnv_profile}"
	mv "${tumor_normal_sample_id}/plot.cnv.png" "${accucopy_cnv_png}"
	mv "${tumor_normal_sample_id}/plot.tre.png" "${accucopy_tre_png}"
	mv "${tumor_normal_sample_id}/het_snp.tsv.gz" "${accucopy_het_snps}"
	mv "${tumor_normal_sample_id}/all_segments.tsv.gz" "${accucopy_segments}"
  	"""
}

// Accucopy Consensus CNV Prep ~ extract and prepare CNV output for consensus
process consensusCnvPrep_accucopy {
  	tag "${tumor_normal_sample_id}"

  	input:
  	tuple val(tumor_normal_sample_id), path(accucopy_cnv_profile) from accucopy_cnv_profile_forConsensusPrep

  	output:
  	tuple val(tumor_normal_sample_id), path(accucopy_somatic_cnv_bed), path(accucopy_somatic_alleles_bed) into final_accucopy_cnv_profile_forConsensus
  	tuple val(tumor_normal_sample_id), path(accucopy_subclones_file) into accucopy_subclones_forConsensusSubclones

  	when:
  	params.accucopy == "on"

  	script:
  	accucopy_somatic_cnv_bed = "${tumor_normal_sample_id}.accucopy.somatic.cnv.bed"
  	accucopy_somatic_alleles_bed = "${tumor_normal_sample_id}.accucopy.somatic.alleles.bed"
  	accucopy_subclones_file = "${tumor_normal_sample_id}.accucopy.subclones.txt"
  	"""
  	accucopy_cnv_profile_postprocessor.sh \
  	"${accucopy_cnv_profile}" \
  	"${tumor_normal_sample_id}" \
  	"${accucopy_subclones_file}" \
  	"${accucopy_somatic_cnv_bed}" \
  	"${accucopy_somatic_alleles_bed}"
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
	tuple val(tumor_normal_sample_id), val(tumor_id), val(normal_id), path(manta_somatic_sv_vcf), path(manta_somatic_sv_vcf_index) into manta_sv_vcf_forPostprocessing
	tuple path(unfiltered_sv_vcf), path(unfiltered_sv_vcf_index)
	tuple path(germline_sv_vcf), path(germline_sv_vcf_index)

	when:
	params.manta == "on"

	script:
	tumor_id = "${tumor_bam.baseName}".replaceFirst(/\..*$/, "")
	normal_id = "${normal_bam.baseName}".replaceFirst(/\..*$/, "")
	tumor_normal_sample_id = "${tumor_id}_vs_${normal_id}"
	zipped_wgs_bed_forManta = "${gatk_bundle_wgs_bed_forManta}.gz"
	unfiltered_sv_vcf = "${tumor_normal_sample_id}.manta.somatic.sv.unfiltered.vcf.gz"
	manta_somatic_config = "${tumor_normal_sample_id}.manta.somatic.config.ini"
	unfiltered_sv_vcf_index = "${unfiltered_sv_vcf}.tbi"
	germline_sv_vcf = "${tumor_normal_sample_id}.manta.germline.sv.vcf.gz"
	germline_sv_vcf_index = "${germline_sv_vcf}.tbi"
	candidate_indel_vcf = "${tumor_normal_sample_id}.manta.somatic.indels.unfiltered.vcf.gz"
	candidate_indel_vcf_index = "${candidate_indel_vcf}.tbi"
	manta_somatic_sv_vcf = "${tumor_normal_sample_id}.manta.somatic.sv.unprocessed.vcf.gz"
	manta_somatic_sv_vcf_index = "${manta_somatic_sv_vcf}.tbi"
	"""
	bgzip < "${gatk_bundle_wgs_bed_forManta}" > "${zipped_wgs_bed_forManta}"
	tabix "${zipped_wgs_bed_forManta}"

	touch "${manta_somatic_config}"
	cat \${MANTA_DIR}/bin/configManta.py.ini \
	| \
	sed 's|enableRemoteReadRetrievalForInsertionsInCancerCallingModes = 0|enableRemoteReadRetrievalForInsertionsInCancerCallingModes = 1|' \
	| \
	sed 's|minPassSomaticScore = 30|minPassSomaticScore = 35|' \
	| \
	sed 's|minCandidateSpanningCount = 3|minCandidateSpanningCount = 4|' >> "${manta_somatic_config}"

	python \${MANTA_DIR}/bin/configManta.py \
	--tumorBam "${tumor_bam}" \
	--normalBam "${normal_bam}" \
	--referenceFasta "${reference_genome_fasta_forManta}" \
	--callRegions "${zipped_wgs_bed_forManta}" \
	--config "${manta_somatic_config}" \
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
	grep -E "^#|PASS" > somaticSV.passonly.vcf

	\${MANTA_DIR}/libexec/convertInversion.py \
	\${MANTA_DIR}/libexec/samtools \
	"${reference_genome_fasta_forManta}" \
	somaticSV.passonly.vcf \
	| \
	bgzip > "${manta_somatic_sv_vcf}"
	tabix "${manta_somatic_sv_vcf}"
	"""
}

// BCFtools filter / view ~ filter out additional false positives based on alternative variant reads in normal sample, prepare VCF for SURVIVOR
process filterAndPostprocessMantaVcf_bcftools {
    tag "${tumor_normal_sample_id}"

    input:
    tuple val(tumor_normal_sample_id), val(tumor_id), val(normal_id), path(manta_somatic_sv_vcf), path(manta_somatic_sv_vcf_index) from manta_sv_vcf_forPostprocessing

    output:
    tuple val(tumor_normal_sample_id), val(tumor_id), path(final_manta_somatic_sv_vcf) into manta_sv_vcf_forSurvivor
    tuple val(tumor_normal_sample_id), path(final_manta_somatic_sv_read_support) into manta_sv_read_support_forAnnotation

    when:
    params.manta == "on"

    script:
    final_manta_somatic_sv_vcf = "${tumor_normal_sample_id}.manta.somatic.sv.vcf"
    final_manta_somatic_sv_read_support = "${tumor_normal_sample_id}.manta.somatic.sv.readsupp.txt"
    """
    touch name.txt
    echo "${normal_id}" >> name.txt

    bcftools filter \
    --output-type v \
    --exclude 'FORMAT/SR[@name.txt:1]>2 || FORMAT/PR[@name.txt:1]>2' \
    "${manta_somatic_sv_vcf}" \
    | \
    bcftools view \
    --output-type v \
    --samples "${tumor_id}" \
    --output-file "${final_manta_somatic_sv_vcf}"

    bcftools query \
    --format '%ID\t[%PR{1}]\t[%SR{1}]\n' \
    --output "${final_manta_somatic_sv_read_support}" \
    "${final_manta_somatic_sv_vcf}"
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


// ~~~~~~~~~~~~~~~~~ SvABA ~~~~~~~~~~~~~~~~~ \\
// START

// Combine all needed reference FASTA files, WGS BED, and reference VCF files into one channel for use in SvABA process
bwa_ref_genome_files.collect()
	.combine( reference_genome_fasta_dict_forSvaba )
	.combine( gatk_bundle_wgs_bed_blacklist_0based_forSvaba )
	.set{ bwa_ref_genome_and_wgs_bed }

bwa_ref_genome_and_wgs_bed.combine( dbsnp_known_indel_ref_vcf )
	.combine( dbsnp_known_indel_ref_vcf_index )
	.combine( simple_and_centromeric_repeats_bed_forSvaba )
	.set{ bwa_ref_genome_wgs_bed_and_ref_files }

// SvABA ~ detecting structural variants using genome-wide local assembly
process svAndIndelCalling_svaba {
	publishDir "${params.output_dir}/somatic/svaba", mode: 'copy', pattern: '*.{somatic.sv.unprocessed.vcf.gz,somatic.sv.unprocessed.vcf.gz.tbi}'
	tag "${tumor_normal_sample_id}"

	input:
	tuple path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path("Homo_sapiens_assembly38.fasta"), path("Homo_sapiens_assembly38.fasta.fai"), path("Homo_sapiens_assembly38.fasta.64.alt"), path("Homo_sapiens_assembly38.fasta.64.amb"), path("Homo_sapiens_assembly38.fasta.64.ann"), path("Homo_sapiens_assembly38.fasta.64.bwt"), path("Homo_sapiens_assembly38.fasta.64.pac"), path("Homo_sapiens_assembly38.fasta.64.sa"), path(reference_genome_fasta_dict_forSvaba), path(gatk_bundle_wgs_bed_blacklist_0based_forSvaba), path(dbsnp_known_indel_ref_vcf), path(dbsnp_known_indel_ref_vcf_index), path(simple_and_centromeric_repeats_bed_forSvaba) from tumor_normal_pair_forSvaba.combine(bwa_ref_genome_wgs_bed_and_ref_files)

	output:
	tuple val(tumor_normal_sample_id), path(filtered_somatic_indel_vcf), path(filtered_somatic_indel_vcf_index) into filtered_indel_vcf_forSvabaBcftools
	tuple val(tumor_normal_sample_id), val(tumor_id), path(svaba_somatic_sv_vcf), path(svaba_somatic_sv_vcf_index), path(svaba_somatic_sv_unclassified_vcf), path(sample_renaming_file) into svaba_sv_vcf_forPostprocessing
	path unfiltered_somatic_indel_vcf
	path unfiltered_somatic_sv_vcf
	path germline_indel_vcf
	path germline_sv_vcf
	path contig_alignment_plot

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
	svaba_somatic_sv_unclassified_vcf = "${tumor_normal_sample_id}.svaba.somatic.sv.unclassified.vcf"
	svaba_somatic_sv_vcf = "${tumor_normal_sample_id}.svaba.somatic.sv.unprocessed.vcf.gz"
	svaba_somatic_sv_vcf_index = "${svaba_somatic_sv_vcf}.tbi"
	sample_renaming_file = "sample_renaming_file.txt"
	"""
	svaba run \
	-t "${tumor_bam}" \
	-n "${normal_bam}" \
	--reference-genome Homo_sapiens_assembly38.fasta \
	--blacklist "${gatk_bundle_wgs_bed_blacklist_0based_forSvaba}" \
	--id-string "${tumor_normal_sample_id}" \
	--dbsnp-vcf "${dbsnp_known_indel_ref_vcf}" \
	--simple-seq-database "${simple_and_centromeric_repeats_bed_forSvaba}" \
	--threads "${task.cpus}" \
	--verbose 1 \
	--g-zip

	mv "${tumor_normal_sample_id}.alignments.txt.gz" "${contig_alignment_plot}"
	mv "${tumor_normal_sample_id}.svaba.unfiltered.somatic.indel.vcf.gz" "${unfiltered_somatic_indel_vcf}"
	mv "${tumor_normal_sample_id}.svaba.unfiltered.somatic.sv.vcf.gz" "${unfiltered_somatic_sv_vcf}"
	mv "${tumor_normal_sample_id}.svaba.somatic.indel.vcf.gz" "${filtered_somatic_indel_vcf}"

	gunzip -c "${tumor_normal_sample_id}.svaba.somatic.sv.vcf.gz" > "${svaba_somatic_sv_unclassified_vcf}"

	svaba_sv_classifier.py \
	"${svaba_somatic_sv_unclassified_vcf}" \
	| \
	bgzip > "${svaba_somatic_sv_vcf}"

	tabix "${filtered_somatic_indel_vcf}"
	tabix "${svaba_somatic_sv_vcf}"

	touch "${sample_renaming_file}"
	echo "${tumor_bam} ${tumor_id}" >> "${sample_renaming_file}"
	"""
}

// BCFtools filter / reheader / view ~ filter out additional false positives based on overall quality score and support
// read mapping quality, prepare VCF for SURVIVOR 
process filterAndPostprocessSvabaVcf_bcftools {
    tag "${tumor_normal_sample_id}"

    input:
    tuple val(tumor_normal_sample_id), val(tumor_id), path(svaba_somatic_sv_vcf), path(svaba_somatic_sv_vcf_index), path(svaba_somatic_sv_unclassified_vcf), path(sample_renaming_file) from svaba_sv_vcf_forPostprocessing

    output:
    tuple val(tumor_normal_sample_id), val(tumor_id), path(final_svaba_somatic_sv_vcf) into svaba_sv_vcf_forSurvivor
    tuple val(tumor_normal_sample_id), path(final_svaba_somatic_sv_read_support) into svaba_sv_read_support_forAnnotation

    when:
    params.svaba == "on"

    script:
    final_svaba_somatic_sv_vcf = "${tumor_normal_sample_id}.svaba.somatic.sv.vcf"
    final_svaba_somatic_sv_read_support = "${tumor_normal_sample_id}.svaba.somatic.sv.readsupp.txt"
    """
    bcftools filter \
    --output-type v \
    --exclude 'QUAL<6' \
    "${svaba_somatic_sv_vcf}" \
    | \
    bcftools filter \
    --output-type v \
    --include 'INFO/MAPQ=60 || INFO/DISC_MAPQ=60' \
    | \
    bcftools reheader \
    --samples "${sample_renaming_file}" \
    | \
    bcftools view \
    --output-type v \
    --samples "${tumor_id}" \
    --output-file "${final_svaba_somatic_sv_vcf}"

    svaba_interchromosomal_mate_finder.sh \
    "${final_svaba_somatic_sv_vcf}" \
    "${svaba_somatic_sv_unclassified_vcf}" > "${tumor_normal_sample_id}.svaba.missingmates.txt"

    cat "${tumor_normal_sample_id}.svaba.missingmates.txt" >> "${final_svaba_somatic_sv_vcf}"

    bcftools query \
    --format '%ID\t[%DR]\t[%SR]\n' \
    --output "${final_svaba_somatic_sv_read_support}" \
    "${final_svaba_somatic_sv_vcf}"
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
	tuple val(tumor_normal_sample_id), path(filtered_somatic_indel_vcf), path(filtered_somatic_indel_vcf_index), path(reference_genome_fasta_forSvabaBcftools), path(reference_genome_fasta_index_forSvabaBcftools), path(reference_genome_fasta_dict_forSvabaBcftools) from filtered_indel_vcf_forSvabaBcftools.combine(reference_genome_bundle_forSvabaBcftools)

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
	tuple val(tumor_normal_sample_id), val(tumor_id), path(delly_somatic_sv_vcf), path(delly_somatic_sv_vcf_index) into delly_sv_vcf_forPostprocessing
	path delly_germline_sv_vcf
	path delly_germline_sv_vcf_index

	when:
	params.delly == "on"

	script:
	tumor_id = "${tumor_bam.baseName}".replaceFirst(/\..*$/, "")
	normal_id = "${normal_bam.baseName}".replaceFirst(/\..*$/, "")
	tumor_normal_sample_id = "${tumor_id}_vs_${normal_id}"
	delly_somatic_sv_vcf = "${tumor_normal_sample_id}.delly.somatic.sv.unprocessed.vcf.gz"
	delly_somatic_sv_vcf_index = "${delly_somatic_sv_vcf}.tbi"
	delly_germline_sv_vcf = "${tumor_normal_sample_id}.delly.germline.sv.vcf.gz"
	delly_germline_sv_vcf_index = "${delly_germline_sv_vcf}.tbi"
	"""
	delly call \
	--genome "${reference_genome_fasta_forDelly}" \
	--exclude "${gatk_bundle_wgs_bed_blacklist_0based_forDelly}" \
	--outfile "${tumor_normal_sample_id}.delly.sv.unfiltered.bcf" \
	"${tumor_bam}" "${normal_bam}"

	touch samples.tsv
	echo "${tumor_id}\ttumor" >> samples.tsv
	echo "${normal_id}\tcontrol" >> samples.tsv

	delly filter \
	--filter somatic \
	--pass \
	--altaf 0.1 \
	--minsize 51 \
	--coverage 10 \
	--samples samples.tsv \
	--outfile "${tumor_normal_sample_id}.delly.somatic.sv.unprocessed.bcf" \
	"${tumor_normal_sample_id}.delly.sv.unfiltered.bcf"

	bcftools isec \
	--nfiles 1 \
	--complement \
	--write 1 \
	--output-type v \
	"${tumor_normal_sample_id}.delly.sv.unfiltered.bcf" \
	"${tumor_normal_sample_id}.delly.somatic.sv.unprocessed.bcf" \
	| \
	bcftools filter \
	--output-type v \
	--include 'FILTER="PASS"' \
	| \
	bgzip > "${delly_germline_sv_vcf}"

	tabix "${delly_germline_sv_vcf}"

	bcftools view \
	--output-type z \
	"${tumor_normal_sample_id}.delly.somatic.sv.unprocessed.bcf" > "${delly_somatic_sv_vcf}"

	tabix "${delly_somatic_sv_vcf}"
	"""
}

// BCFtools filter / reheader / view ~ filter out additional false positives based on overall quality
// score and support read mapping quality, prepare VCF for SURVIVOR 
process filterAndPostprocessDellyVcf_bcftools {
    tag "${tumor_normal_sample_id}"

    input:
    tuple val(tumor_normal_sample_id), val(tumor_id), path(delly_somatic_sv_vcf), path(delly_somatic_sv_vcf_index) from delly_sv_vcf_forPostprocessing

    output:
    tuple val(tumor_normal_sample_id), val(tumor_id), path(final_delly_somatic_sv_vcf) into delly_sv_vcf_forSurvivor
    tuple val(tumor_normal_sample_id), path(final_delly_somatic_sv_read_support) into delly_sv_read_support_forAnnotation

    when:
    params.delly == "on"

    script:
    final_delly_somatic_sv_vcf = "${tumor_normal_sample_id}.delly.somatic.sv.vcf"
    final_delly_somatic_sv_read_support = "${tumor_normal_sample_id}.delly.somatic.sv.readsupp.txt"
    """
    bcftools filter \
    --output-type v \
    --include 'INFO/MAPQ=60 || INFO/SRMAPQ=60' \
    "${delly_somatic_sv_vcf}" \
    | \
    bcftools filter \
    --output-type v \
    --include 'INFO/PE>3 || INFO/SR>3' \
    | \
    bcftools view \
    --output-type v \
    --samples "${tumor_id}" \
    --output-file "${final_delly_somatic_sv_vcf}"

    touch "${tumor_normal_sample_id}.delly.splitmates.txt"
    delly_interchromosomal_record_splitter.sh \
    "${final_delly_somatic_sv_vcf}" > "${tumor_normal_sample_id}.delly.splitmates.txt"

    cat "${tumor_normal_sample_id}.delly.splitmates.txt" >> "${final_delly_somatic_sv_vcf}"

    bcftools query \
    --format '%ID\t%PE\t%SR\n' \
    --output "${final_delly_somatic_sv_read_support}" \
    "${final_delly_somatic_sv_vcf}"
    """
}

// END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \\


// ~~~~~~~~~~~~~~~~ IgCaller ~~~~~~~~~~~~~~~ \\
// START

// Combine all reference FASTA files into one channel for use in IgCaller process
reference_genome_fasta_forIgCaller.combine( reference_genome_fasta_index_forIgCaller )
	.combine( reference_genome_fasta_dict_forIgCaller )
	.set{ reference_genome_bundle_forIgCaller }

// IgCaller ~ characterize the immunoglobulin gene rearrangements and oncogenic translocations
process igRearrangementsAndTranslocations_igcaller {
    publishDir "${params.output_dir}/somatic/igcaller", mode: 'copy', pattern: '*.{tsv}'
    tag "${tumor_normal_sample_id}"

    input:
    tuple path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(reference_genome_fasta_forIgCaller), path(reference_genome_fasta_index_forIgCaller), path(reference_genome_fasta_dict_forIgCaller) from tumor_normal_pair_forIgCaller.combine(reference_genome_bundle_forIgCaller)
     
    output:
    path igcaller_csr_tsv
    path igcaller_igk_tsv
    path igcaller_igl_tsv
    path igcaller_igh_tsv
    path igcaller_filtered_calls_tsv
    path igcaller_oncogenic_rearrangements_tsv

    when:
    params.igcaller == "on"

    script:
    tumor_id = "${tumor_bam.baseName}".replaceFirst(/\..*$/, "")
	normal_id = "${normal_bam.baseName}".replaceFirst(/\..*$/, "")
	tumor_normal_sample_id = "${tumor_id}_vs_${normal_id}"
	igcaller_csr_tsv = "${tumor_normal_sample_id}.igcaller.csr.tsv"
	igcaller_igk_tsv = "${tumor_normal_sample_id}.igcaller.igk.tsv"
	igcaller_igl_tsv = "${tumor_normal_sample_id}.igcaller.igl.tsv"
	igcaller_igh_tsv = "${tumor_normal_sample_id}.igcaller.igh.tsv"
	igcaller_filtered_calls_tsv = "${tumor_normal_sample_id}.igcaller.filtered.tsv"
	igcaller_oncogenic_rearrangements_tsv = "${tumor_normal_sample_id}.igcaller.oncogenic.rearrangements.tsv"
    """
    python3 \${IGCALLER_DIR}/IgCaller.py \
    --inputsFolder \${IGCALLER_DIR}/IgCaller_reference_files/ \
    --genomeVersion hg38 \
    --chromosomeAnnotation ucsc \
    --bamN "${normal_bam}" \
    --bamT "${tumor_bam}" \
    --refGenome "${reference_genome_fasta_forIgCaller}" \
    --outputPath . \
    -@ ${task.cpus}

    mv "${tumor_bam.baseName}_IgCaller/${tumor_bam.baseName}_output_CSR.tsv" "${igcaller_csr_tsv}"
    mv "${tumor_bam.baseName}_IgCaller/${tumor_bam.baseName}_output_IGK.tsv" "${igcaller_igk_tsv}"
    mv "${tumor_bam.baseName}_IgCaller/${tumor_bam.baseName}_output_IGL.tsv" "${igcaller_igl_tsv}"
    mv "${tumor_bam.baseName}_IgCaller/${tumor_bam.baseName}_output_IGH.tsv" "${igcaller_igh_tsv}"
    mv "${tumor_bam.baseName}_IgCaller/${tumor_bam.baseName}_output_filtered.tsv" "${igcaller_filtered_calls_tsv}"
    mv "${tumor_bam.baseName}_IgCaller/${tumor_bam.baseName}_output_oncogenic_IG_rearrangements.tsv" "${igcaller_oncogenic_rearrangements_tsv}"
    """
}

// END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \\


// ~~~~~~~~~ CONSENSUS SNV/INDEL ~~~~~~~~~~~ \\
// START

// MergeVCF ~ merge VCF files by calls, labelling calls by the callers that made them to generate consensus
process mergeAndGenerateConsensusSnvCalls_mergevcf {
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(final_varscan_snv_vcf), path(final_varscan_snv_vcf_index), path(final_mutect_snv_vcf), path(final_mutect_snv_vcf_index), path(final_strelka_snv_vcf), path(final_strelka_snv_vcf_index) from final_varscan_snv_vcf_forConsensus.join(final_mutect_snv_vcf_forConsensus).join(final_strelka_snv_vcf_forConsensus)

	output:
	tuple val(tumor_normal_sample_id), path(consensus_somatic_snv_nosamples_badheader_noformat_vcf), path(snv_vcf_base_header) into consensus_snv_vcf_forConsensusSnvMpileup

	when:
	params.varscan == "on" && params.mutect == "on" && params.strelka == "on"

	script:
	consensus_somatic_snv_nosamples_badheader_noformat_vcf = "${tumor_normal_sample_id}.consensus.somatic.snv.nosamples.badheader.noformat.vcf"
	snv_vcf_base_header = "snv_vcf_base_header.txt"
	"""
	mergevcf \
	--labels varscan,mutect,strelka \
	--ncallers \
	--mincallers 2 \
	"${final_varscan_snv_vcf}" \
	"${final_mutect_snv_vcf}" \
	"${final_strelka_snv_vcf}" \
	| \
	awk '\$1 ~ /^#/ {print \$0;next} {print \$0 | "sort -k1,1V -k2,2n"}' > "${consensus_somatic_snv_nosamples_badheader_noformat_vcf}"

	touch "${snv_vcf_base_header}"
	grep '##fileformat=' "${consensus_somatic_snv_nosamples_badheader_noformat_vcf}" >> "${snv_vcf_base_header}"
	zgrep '##fileDate=' "${final_strelka_snv_vcf}" >> "${snv_vcf_base_header}"
	echo '##source=varscan,mutect,strelka' >> "${snv_vcf_base_header}"
	zgrep '##normal_sample=' "${final_mutect_snv_vcf}" >> "${snv_vcf_base_header}"
	zgrep '##tumor_sample=' "${final_mutect_snv_vcf}" >> "${snv_vcf_base_header}"
	zgrep '##reference=' "${final_strelka_snv_vcf}" | sed -E 's|file.*/||' >> "${snv_vcf_base_header}"
	zgrep '##contig=<ID=chr' "${final_mutect_snv_vcf}" >> "${snv_vcf_base_header}"
	"""
}

// MergeVCF ~ merge VCF files by calls, labelling calls by the callers that made them to generate consensus
process mergeAndGenerateConsensusIndelCalls_mergevcf {
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
	params.varscan == "on" && params.mutect == "on" && params.strelka == "on"

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
	params.varscan == "on" && params.mutect == "on" && params.strelka == "on"

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
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), val(tumor_id), val(normal_id), path(mpileup_supported_consensus_somatic_snv_noformat_vcf), path(snv_mpileup_info_dp_metrics), path(snv_mpileup_normal_format_metrics), path(snv_mpileup_normal_format_metrics_index), path(snv_mpileup_tumor_format_metrics), path(snv_mpileup_tumor_format_metrics_index) from consensus_snv_vcf_forAddFormat.join(consensus_snv_mpileup_metrics_forAddFormat)

	output:
	tuple val(tumor_normal_sample_id), path(snv_consensus_vcf), path(snv_consensus_vcf_index), path(snv_strand_metrics) into consensus_snv_forBedFilters

	when:
	params.varscan == "on" && params.mutect == "on" && params.strelka == "on"

	script:
	snv_consensus_vcf_info_header = "snv_consensus_vcf_info_header.txt"
	snv_consensus_vcf_format_headers = "snv_consensus_vcf_format_header.txt"
	snv_consensus_vcf = "${tumor_normal_sample_id}.consensus.somatic.snv.vcf.gz"
	snv_consensus_vcf_index = "${snv_consensus_vcf}.tbi"
	snv_strand_metrics = "${tumor_normal_sample_id}.snv.strand.metrics.txt"
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

	touch "${snv_consensus_vcf_info_header}"
	echo '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth across samples (normal sample DP + tumor sample DP)">' >> "${snv_consensus_vcf_info_header}"
	echo '##INFO=<ID=AC,Number=1,Type=Integer,Description="Count of ALT allele reads in tumor sample">' >> "${snv_consensus_vcf_info_header}"
	echo '##INFO=<ID=VAF,Number=1,Type=Float,Description="Variant allele frequency, expressed as fraction of ALT allele reads in total read depth in tumor sample (tumor sample ALT AC / tumor sample DP)">' >> "${snv_consensus_vcf_info_header}"

	bcftools annotate \
	--output-type z \
	--annotations "${tumor_normal_sample_id}.snv.mpileup.info.metrics.txt.gz" \
	--header-lines "${snv_consensus_vcf_info_header}" \
	--columns CHROM,POS,REF,ALT,INFO/DP,INFO/AC,INFO/VAF \
	--output "${tumor_normal_sample_id}.ms.consensus.somatic.snv.info.noformat.vcf.gz" \
	"${mpileup_supported_consensus_somatic_snv_noformat_vcf}"

	touch "${snv_consensus_vcf_format_headers}"
	echo '##FORMAT=<ID=DPS,Number=1,Type=Integer,Description="Total read depth in sample">' >> "${snv_consensus_vcf_format_headers}"
	echo '##FORMAT=<ID=ACS,Number=R,Type=Integer,Description="Count of REF,ALT allele reads in sample">' >> "${snv_consensus_vcf_format_headers}"
	echo '##FORMAT=<ID=ACFS,Number=R,Type=Integer,Description="Count of REF,ALT allele reads on forward(+) strand in sample">' >> "${snv_consensus_vcf_format_headers}"
	echo '##FORMAT=<ID=ACRS,Number=R,Type=Integer,Description="Count of REF,ALT allele reads on reverse(-) strand in sample">' >> "${snv_consensus_vcf_format_headers}"

	bcftools annotate \
	--output-type z \
	--samples "${normal_id}" \
	--annotations "${snv_mpileup_normal_format_metrics}" \
	--header-lines "${snv_consensus_vcf_format_headers}" \
	--columns CHROM,POS,REF,ALT,FORMAT/DPS,FORMAT/ACS,FORMAT/ACFS,FORMAT/ACRS \
	--output "${tumor_normal_sample_id}.ms.consensus.somatic.snv.info.halfformat.vcf.gz" \
	"${tumor_normal_sample_id}.ms.consensus.somatic.snv.info.noformat.vcf.gz"

	bcftools annotate \
	--output-type z \
	--samples "${tumor_id}" \
	--annotations "${snv_mpileup_tumor_format_metrics}" \
	--header-lines "${snv_consensus_vcf_format_headers}" \
	--columns CHROM,POS,REF,ALT,FORMAT/DPS,FORMAT/ACS,FORMAT/ACFS,FORMAT/ACRS \
	--remove FORMAT/GT \
	--output "${tumor_normal_sample_id}.ms.consensus.somatic.snv.info.format.vcf.gz" \
	"${tumor_normal_sample_id}.ms.consensus.somatic.snv.info.halfformat.vcf.gz"

	bcftools filter \
	--output-type v \
	--exclude 'INFO/AC<3 | INFO/VAF<0.01' \
	"${tumor_normal_sample_id}.ms.consensus.somatic.snv.info.format.vcf.gz" \
	| \
	bcftools filter \
	--output-type v \
	--exclude 'FILTER="LOWSUPPORT"' - \
	| \
	grep -v '##bcftools_annotate' \
	| \
	bgzip > "${snv_consensus_vcf}"
	tabix "${snv_consensus_vcf}"

	bcftools query \
	--format '%CHROM\t%POS\t[%ACFS]\t[%ACRS]\n' \
	--samples "${tumor_id}" \
	--output "${snv_strand_metrics}" \
	"${snv_consensus_vcf}"
	"""
}

// VCFtools ~ final hard filters of SNVs based on simple/centromeric repeats and strand bias BED files before annotation of SNVs
process repeatsAndStrandBiasFilterSnvs_vcftools {
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(snv_consensus_vcf), path(snv_consensus_vcf_index), path(snv_strand_metrics), path(simple_and_centromeric_repeats_bed) from consensus_snv_forBedFilters.combine(simple_and_centromeric_repeats_bed_forSnvBedFilter)

	output:
	tuple val(tumor_normal_sample_id), path(hq_snv_consensus_vcf), path(hq_snv_consensus_vcf_index) into high_quality_consensus_snv_forAnnotation

	when:
	params.varscan == "on" && params.mutect == "on" && params.strelka == "on"

	script:
	strand_bias_filter_bed = "${tumor_normal_sample_id}.snv.strandbias.bed"
	hq_snv_consensus_vcf = "${tumor_normal_sample_id}.hq.consensus.somatic.snv.vcf.gz"
	hq_snv_consensus_vcf_index = "${hq_snv_consensus_vcf}.tbi"
	"""
	Rscript --vanilla  \
	${workflow.projectDir}/bin/strand_bias_proportion_tester.R \
	"${snv_strand_metrics}" \
	"${strand_bias_filter_bed}"

	vcftools \
	--gzvcf "${snv_consensus_vcf}" \
	--exclude-bed "${strand_bias_filter_bed}" \
	--recode \
	--recode-INFO-all \
	--stdout \
	| \
	vcftools \
	--vcf - \
	--exclude-bed "${simple_and_centromeric_repeats_bed}" \
	--recode \
	--recode-INFO-all \
	--stdout \
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
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), val(tumor_id), val(normal_id), path(mpileup_supported_consensus_somatic_indel_noformat_vcf), path(indel_mpileup_info_dp_metrics), path(indel_mpileup_normal_format_metrics), path(indel_mpileup_normal_format_metrics_index), path(indel_mpileup_tumor_format_metrics), path(indel_mpileup_tumor_format_metrics_index) from consensus_indel_vcf_forAddFormat.join(consensus_indel_mpileup_metrics_forAddFormat)

	output:
	tuple val(tumor_normal_sample_id), path(indel_consensus_vcf), path(indel_consensus_vcf_index), path(indel_strand_metrics) into consensus_indel_forBedFilters

	when:
	params.varscan == "on" && params.mutect == "on" && params.strelka == "on" && params.svaba == "on"

	script:
	indel_consensus_vcf_info_header = "indel_consensus_vcf_info_header.txt"
	indel_consensus_vcf_format_headers = "indel_consensus_vcf_format_header.txt"
	indel_consensus_vcf = "${tumor_normal_sample_id}.consensus.somatic.indel.vcf.gz"
	indel_consensus_vcf_index = "${indel_consensus_vcf}.tbi"
	indel_strand_metrics = "${tumor_normal_sample_id}.indel.strand.mertrics.txt"
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

	touch "${indel_consensus_vcf_info_header}"
	echo '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth across samples (normal sample DP + tumor sample DP)">' >> "${indel_consensus_vcf_info_header}"
	echo '##INFO=<ID=AC,Number=1,Type=Integer,Description="Count of ALT allele reads in tumor sample">' >> "${indel_consensus_vcf_info_header}"
	echo '##INFO=<ID=VAF,Number=1,Type=Float,Description="Variant allele frequency, expressed as fraction of ALT allele reads in total read depth in tumor sample (tumor sample ALT AC / tumor sample DP)">' >> "${indel_consensus_vcf_info_header}"

	bcftools annotate \
	--output-type z \
	--annotations "${tumor_normal_sample_id}.indel.mpileup.info.metrics.txt.gz" \
	--header-lines "${indel_consensus_vcf_info_header}" \
	--columns CHROM,POS,REF,ALT,INFO/DP,INFO/AC,INFO/VAF \
	--output "${tumor_normal_sample_id}.ms.consensus.somatic.indel.info.noformat.vcf.gz" \
	"${mpileup_supported_consensus_somatic_indel_noformat_vcf}"

	touch "${indel_consensus_vcf_format_headers}"
	echo '##FORMAT=<ID=DPS,Number=1,Type=Integer,Description="Total read depth in sample">' >> "${indel_consensus_vcf_format_headers}"
	echo '##FORMAT=<ID=ACS,Number=R,Type=Integer,Description="Count of REF,ALT allele reads in sample">' >> "${indel_consensus_vcf_format_headers}"
	echo '##FORMAT=<ID=ACFS,Number=R,Type=Integer,Description="Count of REF,ALT allele reads on forward(+) strand in sample">' >> "${indel_consensus_vcf_format_headers}"
	echo '##FORMAT=<ID=ACRS,Number=R,Type=Integer,Description="Count of REF,ALT allele reads on reverse(-) strand in sample">' >> "${indel_consensus_vcf_format_headers}"

	bcftools annotate \
	--output-type z \
	--samples "${normal_id}" \
	--annotations "${indel_mpileup_normal_format_metrics}" \
	--header-lines "${indel_consensus_vcf_format_headers}" \
	--columns CHROM,POS,REF,ALT,FORMAT/DPS,FORMAT/ACS,FORMAT/ACFS,FORMAT/ACRS \
	--output "${tumor_normal_sample_id}.ms.consensus.somatic.indel.info.halfformat.vcf.gz" \
	"${tumor_normal_sample_id}.ms.consensus.somatic.indel.info.noformat.vcf.gz"

	bcftools annotate \
	--output-type z \
	--samples "${tumor_id}" \
	--annotations "${indel_mpileup_tumor_format_metrics}" \
	--header-lines "${indel_consensus_vcf_format_headers}" \
	--columns CHROM,POS,REF,ALT,FORMAT/DPS,FORMAT/ACS,FORMAT/ACFS,FORMAT/ACRS \
	--remove FORMAT/GT \
	--output "${tumor_normal_sample_id}.ms.consensus.somatic.indel.info.format.vcf.gz" \
	"${tumor_normal_sample_id}.ms.consensus.somatic.indel.info.halfformat.vcf.gz"

	bcftools filter \
	--output-type v \
	--exclude 'INFO/AC<3 | INFO/VAF<0.01' \
	"${tumor_normal_sample_id}.ms.consensus.somatic.indel.info.format.vcf.gz" \
	| \
	bcftools filter \
	--output-type v \
	--exclude 'FILTER="LOWSUPPORT"' - \
	| \
	grep -v '##bcftools_annotate' \
	| \
	bgzip > "${indel_consensus_vcf}"
	tabix "${indel_consensus_vcf}"

	bcftools query \
	--format '%CHROM\t%POS\t[%ACFS]\t[%ACRS]\n' \
	--samples "${tumor_id}" \
	--output "${indel_strand_metrics}" \
	"${indel_consensus_vcf}"
	"""
}

// VCFtools ~ final hard filters of InDels based on simple/centromeric repeats and strand bias BED files before annotation of SNVs
process repeatsAndStrandBiasFilterIndels_vcftools {
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(indel_consensus_vcf), path(indel_consensus_vcf_index), path(indel_strand_metrics), path(simple_and_centromeric_repeats_bed) from consensus_indel_forBedFilters.combine(simple_and_centromeric_repeats_bed_forIndelBedFilter)

	output:
	tuple val(tumor_normal_sample_id), path(hq_indel_consensus_vcf), path(hq_indel_consensus_vcf_index) into high_quality_consensus_indel_forAnnotation

	when:
	params.varscan == "on" && params.mutect == "on" && params.strelka == "on" && params.svaba == "on"

	script:
	strand_bias_filter_bed = "${tumor_normal_sample_id}.indel.strandbias.bed"
	hq_indel_consensus_vcf = "${tumor_normal_sample_id}.hq.consensus.somatic.indel.vcf.gz"
	hq_indel_consensus_vcf_index = "${hq_indel_consensus_vcf}.tbi"
	"""
	Rscript --vanilla  \
	${workflow.projectDir}/bin/strand_bias_proportion_tester.R \
	"${indel_strand_metrics}" \
	"${strand_bias_filter_bed}"

	vcftools \
	--gzvcf "${indel_consensus_vcf}" \
	--exclude-bed "${strand_bias_filter_bed}" \
	--recode \
	--recode-INFO-all \
	--stdout \
	| \
	vcftools \
	--vcf - \
	--exclude-bed "${simple_and_centromeric_repeats_bed}" \
	--recode \
	--recode-INFO-all \
	--stdout \
	| \
	bgzip > "${hq_indel_consensus_vcf}"

	tabix "${hq_indel_consensus_vcf}"
	"""
}

// END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \\


// ~~~~~~~~~~ CONSENSUS CNV BED ~~~~~~~~~~~~ \\
// START

// Determine which CNV consensus mechanism to use based on number of tools included, either 3 or 4 
if( params.battenberg == "on" & params.controlfreec == "on" & params.sclust == "on" & params.accucopy == "on" ) {
  	cnv_output_forFourWayConsensus = final_battenberg_cnv_profile_forConsensus.join(final_controlfreec_cnv_profile_forConsensus).join(final_sclust_cnv_profile_forConsensus).join(final_accucopy_cnv_profile_forConsensus)
  	cnv_output_forThreeWayConsensus = Channel.empty()

} else if( params.battenberg == "on" & params.controlfreec == "on" & params.sclust == "off" & params.accucopy == "on" ) {
  	cnv_output_forFourWayConsensus = Channel.empty()
  	cnv_output_forThreeWayConsensus = final_battenberg_cnv_profile_forConsensus.join(final_controlfreec_cnv_profile_forConsensus).join(final_accucopy_cnv_profile_forConsensus)

} else {
  	cnv_output_forFourWayConsensus = Channel.empty()
  	cnv_output_forThreeWayConsensus = Channel.empty()
}

// BEDtools unionbedg 4-way ~ transform CNV output into BED files then generate merged CNV segment file
process fourWayMergeAndGenerateConsensusCnvCalls_bedtools {
  	tag "${tumor_normal_sample_id}"

  	input:
  	tuple val(tumor_normal_sample_id), path(battenberg_somatic_cnv_bed), path(battenberg_somatic_alleles_bed), path(control_freec_somatic_cnv_bed), path(control_freec_somatic_alleles_bed), path(sclust_somatic_cnv_bed), path(sclust_somatic_alleles_bed), path(accucopy_somatic_cnv_bed), path(accucopy_somatic_alleles_bed) from cnv_output_forFourWayConsensus

  	output:
  	tuple val(tumor_normal_sample_id), val(four_way_consensus_mechanism), path(four_way_consensus_merged_cnv_alleles_bed) into four_way_consensus_cnv_and_allele_bed_forConsensusCnvTransform

  	when:
  	params.battenberg == "on" & params.controlfreec == "on" & params.sclust == "on" & params.accucopy == "on"

  	script:
  	four_way_consensus_mechanism = "four_way"
  	four_way_merged_cnv_bed = "${tumor_normal_sample_id}.merged.somatic.cnv.bed"
  	four_way_merged_alleles_bed = "${tumor_normal_sample_id}.merged.somatic.alleles.bed"
  	four_way_consensus_cnv_bed = "${tumor_normal_sample_id}.consensus.somatic.cnv.bed"
  	four_way_consensus_alleles_bed = "${tumor_normal_sample_id}.consensus.somatic.alleles.bed"
  	four_way_consensus_merged_cnv_alleles_bed = "${tumor_normal_sample_id}.consensus.somatic.cnv.alleles.merged.bed"
  	"""
  	### Create consensus total copy number file ###
  	bedtools unionbedg \
  	-filler . \
  	-i "${battenberg_somatic_cnv_bed}" "${control_freec_somatic_cnv_bed}" "${sclust_somatic_cnv_bed}" "${accucopy_somatic_cnv_bed}" \
  	-header \
  	-names battenberg_total_cn controlfreec_total_cn sclust_total_cn accucopy_total_cn > "${four_way_merged_cnv_bed}"

  	four_way_consensus_cnv_generator.py \
  	<(grep -v 'chrom' "${four_way_merged_cnv_bed}") \
  	"${four_way_consensus_cnv_bed}"

  	### Create consensus major and minor allele file ###
  	bedtools unionbedg \
  	-filler . \
  	-i "${battenberg_somatic_alleles_bed}" "${control_freec_somatic_alleles_bed}" "${sclust_somatic_alleles_bed}" "${accucopy_somatic_alleles_bed}" \
  	-header \
  	-names battenberg_major_minor_alleles controlfreec_major_minor_alleles sclust_major_minor_alleles accucopy_major_minor_alleles > "${four_way_merged_alleles_bed}"

  	four_way_consensus_allele_generator.py \
  	<(grep -v 'chrom' "${four_way_merged_alleles_bed}") \
  	"${four_way_consensus_alleles_bed}"

  	### Merge both consensus CNV and called alleles per segment ###
  	paste "${four_way_consensus_cnv_bed}" <(cut -f 4-10 "${four_way_consensus_alleles_bed}") \
  	| \
  	awk 'BEGIN {OFS="\t"} {print \$1,\$2,\$3,\$4,\$10,\$11,\$5,\$12,\$6,\$13,\$7,\$14,\$8,\$15,\$9,\$16}' > "${four_way_consensus_merged_cnv_alleles_bed}"
  	"""
}

// BEDtools unionbedg 3-way ~ transform CNV output into BED files then generate merged CNV segment file
process threeWayMergeAndGenerateConsensusCnvCalls_bedtools {
  	tag "${tumor_normal_sample_id}"

  	input:
  	tuple val(tumor_normal_sample_id), path(battenberg_somatic_cnv_bed), path(battenberg_somatic_alleles_bed), path(control_freec_somatic_cnv_bed), path(control_freec_somatic_alleles_bed), path(accucopy_somatic_cnv_bed), path(accucopy_somatic_alleles_bed) from cnv_output_forThreeWayConsensus

  	output:
  	tuple val(tumor_normal_sample_id), val(three_way_consensus_mechanism), path(three_way_consensus_merged_cnv_alleles_bed) into three_way_consensus_cnv_and_allele_bed_forConsensusCnvTransform

  	when:
  	params.battenberg == "on" & params.controlfreec == "on" & params.sclust == "off" & params.accucopy == "on"

  	script:
  	three_way_consensus_mechanism = "three_way"
  	three_way_merged_cnv_bed = "${tumor_normal_sample_id}.merged.somatic.cnv.bed"
  	three_way_merged_alleles_bed = "${tumor_normal_sample_id}.merged.somatic.alleles.bed"
  	three_way_consensus_cnv_bed = "${tumor_normal_sample_id}.consensus.somatic.cnv.bed"
  	three_way_consensus_alleles_bed = "${tumor_normal_sample_id}.consensus.somatic.alleles.bed"
  	three_way_consensus_merged_cnv_alleles_bed = "${tumor_normal_sample_id}.consensus.somatic.cnv.alleles.merged.bed"
  	"""
  	### Create consensus total copy number file ###
  	bedtools unionbedg \
  	-filler . \
  	-i "${battenberg_somatic_cnv_bed}" "${control_freec_somatic_cnv_bed}" "${accucopy_somatic_cnv_bed}" \
  	-header \
  	-names battenberg_total_cn controlfreec_total_cn accucopy_total_cn > "${three_way_merged_cnv_bed}"

  	three_way_consensus_cnv_generator.py \
  	<(grep -v 'chrom' "${three_way_merged_cnv_bed}") \
  	"${three_way_consensus_cnv_bed}"

  	### Create consensus major and minor allele file ###
  	bedtools unionbedg \
  	-filler . \
  	-i "${battenberg_somatic_alleles_bed}" "${control_freec_somatic_alleles_bed}" "${accucopy_somatic_alleles_bed}" \
  	-header \
  	-names battenberg_major_minor_alleles controlfreec_major_minor_alleles accucopy_major_minor_alleles > "${three_way_merged_alleles_bed}"

  	three_way_consensus_allele_generator.py \
  	<(grep -v 'chrom' "${three_way_merged_alleles_bed}") \
  	"${three_way_consensus_alleles_bed}"

  	### Merge both consensus CNV and called alleles per segment ###
  	paste "${three_way_consensus_cnv_bed}" <(cut -f 4-9 "${three_way_consensus_alleles_bed}") \
  	| \
  	awk 'BEGIN {OFS="\t"} {print \$1,\$2,\$3,\$4,\$9,\$10,\$5,\$11,\$6,\$12,\$7,\$13,\$8,\$14}' > "${three_way_consensus_merged_cnv_alleles_bed}"
  	"""
}

// Set the input for the high-quality consensus CNV BED transform process based on 3 or 4 way mechanism
if( params.battenberg == "on" & params.controlfreec == "on" & params.sclust == "on" & params.accucopy == "on" ) {
  	consensus_cnv_and_allele_bed_forConsensusCnvTransform = four_way_consensus_cnv_and_allele_bed_forConsensusCnvTransform

} else if( params.battenberg == "on" & params.controlfreec == "on" & params.sclust == "off" & params.accucopy == "on" ) {
  	consensus_cnv_and_allele_bed_forConsensusCnvTransform = three_way_consensus_cnv_and_allele_bed_forConsensusCnvTransform

} else {
  	consensus_cnv_and_allele_bed_forConsensusCnvTransform = Channel.empty()
}

// Tidyverse ~ transform consensus CNV BED file to concise high quality format
process highQualityTransformConsensusCnvs_tidyverse {
    tag "${tumor_normal_sample_id}"

    input:
    tuple val(tumor_normal_sample_id), val(consensus_mechanism), path(consensus_merged_cnv_alleles_bed), path(sample_sex) from consensus_cnv_and_allele_bed_forConsensusCnvTransform.join(sex_of_sample_forConsensusCnvTransform)

    output:
    tuple val(tumor_normal_sample_id), path(hq_consensus_cnv_bed) into hq_consensus_cnv_bed_forAnnotation

    when:
    params.battenberg == "on" && params.controlfreec == "on" && params.accucopy == "on"

    script:
    hq_consensus_cnv_bed = "${tumor_normal_sample_id}.hq.consensus.somatic.cnv.bed"
    """
    sex=\$(cut -d ' ' -f 1 "${sample_sex}")

    Rscript --vanilla  \
	${workflow.projectDir}/bin/high_quality_cnv_bed_transformer.R \
    "${consensus_merged_cnv_alleles_bed}" \
    \${sex} \
    "${hq_consensus_cnv_bed}" \
    "${consensus_mechanism}"
    """
}

// AnnotSV ~ download the reference files used for AnnotSV annotation, if needed
process downloadAnnotsvAnnotationReferences_annotsv {
    publishDir "references/hg38", mode: 'copy'

    output:
    path cached_ref_dir_annotsv into annotsv_ref_dir_from_process_forSvAnnotation, annotsv_ref_dir_from_process_forCnvAnnotation

    when:
    params.annotsv_ref_cached == "no"

    script:
    cached_ref_dir_annotsv = "annotations_human_annotsv_hg38"
    """
    mkdir -p "${cached_ref_dir_annotsv}" && \
    cd "${cached_ref_dir_annotsv}" && \
    curl -C - -LO https://www.lbgi.fr/~geoffroy/Annotations/Annotations_Human_3.1.1.tar.gz && \
    tar -zxf Annotations_Human_3.1.1.tar.gz && \
    rm Annotations_Human_3.1.1.tar.gz
    """
}

// Depending on whether the reference files used for AnnotSV annotation was pre-downloaded, set the input
// channel for the AnnotSV CNV annotation process
if( params.annotsv_ref_cached == "yes" ) {
    annotsv_ref_dir_forCnvAnnotation = annotsv_ref_dir_pre_downloaded_forCnvAnnotation
}
else {
    annotsv_ref_dir_forCnvAnnotation = annotsv_ref_dir_from_process_forCnvAnnotation
}

// AnnotSV ~ annotate consensus CNV calls with multiple resources
process annotateConsensusCnvCalls_annotsv {
    publishDir "${params.output_dir}/somatic/consensus/${tumor_normal_sample_id}", mode: 'copy', pattern: '*.{bed}'
    tag  "${tumor_normal_sample_id}"

    input:
    tuple val(tumor_normal_sample_id), path(hq_consensus_cnv_bed), path(annotsv_ref_dir_bundle) from hq_consensus_cnv_bed_forAnnotation.combine(annotsv_ref_dir_forCnvAnnotation)

    output:
    path hq_consensus_cnv_annotated_bed

    when:
    params.battenberg == "on" && params.controlfreec == "on" && params.accucopy == "on"

    script:
    hq_consensus_cnv_annotated_bed = "${tumor_normal_sample_id}.hq.consensus.somatic.cnv.annotated.bed"
    """
    \$ANNOTSV/bin/AnnotSV \
    -annotationsDir "${annotsv_ref_dir_bundle}" \
    -annotationMode full \
    -genomeBuild GRCh38 \
    -outputDir . \
    -outputFile "${tumor_normal_sample_id}.hq.consensus.somatic.cnv.annotated" \
    -SVinputFile "${hq_consensus_cnv_bed}" \
    -SVminSize 1 \
    -includeCI 0 \
    -tx ENSEMBL

    touch "${hq_consensus_cnv_annotated_bed}"
    echo "chr\tstart\tend\tconsensus_total_cn\tconsensus_major_allele\tconsensus_minor_allele\ttype\tclass\tallele_status\tconsensus_conf_rating\tcytoband\tgenes" >> "${hq_consensus_cnv_annotated_bed}"
    grep -v 'SV_chrom' "${tumor_normal_sample_id}.hq.consensus.somatic.cnv.annotated.tsv" \
    | \
    awk 'BEGIN {OFS="\t"} {print "chr"\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10,\$11,\$13,\$14}' \
    | \
    sort -k1,1V -k2,2n >> "${hq_consensus_cnv_annotated_bed}"
    """
}

// Determine which subclone consensus mechanism to use based on number of tools included, either 3 or 4 
if( params.battenberg == "on" & params.controlfreec == "on" & params.sclust == "on" & params.accucopy == "on" ) {
    subclone_output_forFourWayConsensus = battenberg_subclones_forConsensusSubclones.join(control_freec_subclones_forConsensusSubclones).join(sclust_subclones_forConsensusSubclones).join(accucopy_subclones_forConsensusSubclones)
    subclone_output_forThreeWayConsensus = Channel.empty()

} else if( params.battenberg == "on" & params.controlfreec == "on" & params.sclust == "off" & params.accucopy == "on" ) {
    subclone_output_forFourWayConsensus = Channel.empty()
    subclone_output_forThreeWayConsensus = battenberg_subclones_forConsensusSubclones.join(control_freec_subclones_forConsensusSubclones).join(accucopy_subclones_forConsensusSubclones)

} else {
    subclone_output_forFourWayConsensus = Channel.empty()
    subclone_output_forThreeWayConsensus = Channel.empty()
}

// 4-way merge four subclone CNV segment calls, simple concatenation of per tool output
process fourWayMergeSubclonalCnvCalls {
	publishDir "${params.output_dir}/somatic/consensus/${tumor_normal_sample_id}", mode: 'copy', pattern: '*.{txt}'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(battenberg_subclones_file), path(control_freec_subclones_file), path(sclust_subclones_file), path(accucopy_subclones_file) from subclone_output_forFourWayConsensus

	output:
	path four_way_consensus_subclonal_cnv_file

	when:
	params.battenberg == "on" && params.controlfreec == "on" && params.sclust == "on" && params.accucopy == "on"

	script:
	four_way_consensus_subclonal_cnv_file = "${tumor_normal_sample_id}.consensus.somatic.cnv.subclonal.txt"
	"""
	touch "${four_way_consensus_subclonal_cnv_file}"

	echo "### Battenberg ###" >> "${four_way_consensus_subclonal_cnv_file}"
	echo "" >> "${four_way_consensus_subclonal_cnv_file}"
	cat "${battenberg_subclones_file}" >> "${four_way_consensus_subclonal_cnv_file}"
	echo "" >> "${four_way_consensus_subclonal_cnv_file}"

	echo "### Control-FREEC ###" >> "${four_way_consensus_subclonal_cnv_file}"
	echo "" >> "${four_way_consensus_subclonal_cnv_file}"
	cat "${control_freec_subclones_file}" >> "${four_way_consensus_subclonal_cnv_file}"
	echo "" >> "${four_way_consensus_subclonal_cnv_file}"

	echo "### Sclust ###" >> "${four_way_consensus_subclonal_cnv_file}"
	echo "" >> "${four_way_consensus_subclonal_cnv_file}"
	cut -f 2-11 "${sclust_subclones_file}" >> "${four_way_consensus_subclonal_cnv_file}"
	echo "" >> "${four_way_consensus_subclonal_cnv_file}"

	echo "### Accucopy ###" >> "${four_way_consensus_subclonal_cnv_file}"
	echo "" >> "${four_way_consensus_subclonal_cnv_file}"
	cat "${accucopy_subclones_file}" >> "${four_way_consensus_subclonal_cnv_file}"
	echo "" >> "${four_way_consensus_subclonal_cnv_file}"
	"""
}

// 3-way merge four subclone CNV segment calls, simple concatenation of per tool output
process threeWayMergeSubclonalCnvCalls {
	publishDir "${params.output_dir}/somatic/consensus/${tumor_normal_sample_id}", mode: 'copy', pattern: '*.{txt}'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(battenberg_subclones_file), path(control_freec_subclones_file), path(accucopy_subclones_file) from subclone_output_forThreeWayConsensus

	output:
	path three_way_consensus_subclonal_cnv_file

	when:
	params.battenberg == "on" && params.controlfreec == "on" && params.sclust == "off" && params.accucopy == "on"

	script:
	three_way_consensus_subclonal_cnv_file = "${tumor_normal_sample_id}.consensus.somatic.cnv.subclonal.txt"
	"""
	touch "${three_way_consensus_subclonal_cnv_file}"

	echo "### Battenberg ###" >> "${three_way_consensus_subclonal_cnv_file}"
	echo "" >> "${three_way_consensus_subclonal_cnv_file}"
	cat "${battenberg_subclones_file}" >> "${three_way_consensus_subclonal_cnv_file}"
	echo "" >> "${three_way_consensus_subclonal_cnv_file}"

	echo "### Control-FREEC ###" >> "${three_way_consensus_subclonal_cnv_file}"
	echo "" >> "${three_way_consensus_subclonal_cnv_file}"
	cat "${control_freec_subclones_file}" >> "${three_way_consensus_subclonal_cnv_file}"
	echo "" >> "${three_way_consensus_subclonal_cnv_file}"

	echo "### Accucopy ###" >> "${three_way_consensus_subclonal_cnv_file}"
	echo "" >> "${three_way_consensus_subclonal_cnv_file}"
	cat "${accucopy_subclones_file}" >> "${three_way_consensus_subclonal_cnv_file}"
	echo "" >> "${three_way_consensus_subclonal_cnv_file}"
	"""
}

// END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \\


// ~~~~~~~~~~~ CONSENSUS SV VCF ~~~~~~~~~~~~ \\
// START

// SURVIVOR ~ merge SV VCF files to generate a consensus
process mergeAndGenerateConsensusSvCalls_survivor {
    tag  "${tumor_normal_sample_id}"

    input:
    tuple val(tumor_normal_sample_id), val(tumor_id), path(final_manta_somatic_sv_vcf), path(final_svaba_somatic_sv_vcf), path(final_delly_somatic_sv_vcf) from manta_sv_vcf_forSurvivor.join(svaba_sv_vcf_forSurvivor, by: [0,1]).join(delly_sv_vcf_forSurvivor, by: [0,1]) 

    output:
    tuple val(tumor_normal_sample_id), val(tumor_id), path(consensus_somatic_sv_badheader_vcf) into consensus_sv_vcf_forFilterPrep

    when:
    params.manta == "on" && params.svaba == "on" && params.delly == "on"

    script:
    consensus_somatic_sv_badheader_vcf = "${tumor_normal_sample_id}.consensus.somatic.sv.badheader.vcf"
    """
    touch input_vcf_list.txt
    ls *.vcf >> input_vcf_list.txt

    SURVIVOR merge \
    input_vcf_list.txt \
    1000 \
    1 \
    0 \
    1 \
    0 \
    51 \
    "${consensus_somatic_sv_badheader_vcf}"
    """
}

// VAtools vcf-genotype-annotator ~ add sample name to SURVIVOR consensus SV VCF for use in Duphold filtering to remove false positives
process prepConsensusSvVcfForFpFiltering_vatools {
    tag "${tumor_normal_sample_id}"

    input:
    tuple val(tumor_normal_sample_id), val(tumor_id), path(consensus_somatic_sv_badheader_vcf) from consensus_sv_vcf_forFilterPrep

    output:
    tuple val(tumor_normal_sample_id), path(consensus_somatic_sv_vcf) into consensus_sv_vcf_forConsensusSvFpFilter

    when:
    params.manta == "on" && params.svaba == "on" && params.delly == "on"

    script:
    consensus_somatic_sv_vcf = "${tumor_normal_sample_id}.consensus.somatic.sv.vcf"
    """
    sed 's|${tumor_id}\t${tumor_id}_1\t${tumor_id}_2|${tumor_id}_delly\t${tumor_id}_manta\t${tumor_id}_svaba|' "${tumor_normal_sample_id}.consensus.somatic.sv.badheader.vcf" \
    | \
    sed 's|SVTYPE=TRA|SVTYPE=BND|' > "${tumor_normal_sample_id}.consensus.somatic.sv.halffixed.vcf"

    vcf-genotype-annotator \
    --output-vcf "${consensus_somatic_sv_vcf}" \
    "${tumor_normal_sample_id}.consensus.somatic.sv.halffixed.vcf" \
    "${tumor_id}" \
    .
    """
}

// Combine all reference FASTA files into one channel for use in duphold process
reference_genome_fasta_forConsensusSvFpFilter.combine( reference_genome_fasta_index_forConsensusSvFpFilter )
    .set{ reference_genome_bundle_forConsensusSvFpFilter }

// duphold ~ efficiently annotate SV calls with sequence depth information to reduce false positive deletion and duplication calls
process falsePostiveSvFiltering_duphold {
    tag "${tumor_normal_sample_id}"

    input:
    tuple val(tumor_normal_sample_id), path(consensus_somatic_sv_vcf), path(tumor_bam), path(tumor_bam_index), path(reference_genome_fasta_forConsensusSvFpFilter), path(reference_genome_fasta_index_forConsensusSvFpFilter) from consensus_sv_vcf_forConsensusSvFpFilter.join(bam_forConsensusSvFpFilter).combine(reference_genome_bundle_forConsensusSvFpFilter)

    output:
    tuple val(tumor_normal_sample_id), path(consensus_somatic_sv_fpmarked_vcf) into consensus_sv_vcf_forFpFiltering

    when:
    params.manta == "on" && params.svaba == "on" && params.delly == "on"

    script:
    consensus_somatic_sv_fpmarked_vcf = "${tumor_normal_sample_id}.consensus.somatic.sv.fpmarked.vcf"
    """
    duphold \
    --vcf "${consensus_somatic_sv_vcf}" \
    --bam "${tumor_bam}" \
    --fasta "${reference_genome_fasta_forConsensusSvFpFilter}" \
    --threads ${task.cpus} \
    --output "${consensus_somatic_sv_fpmarked_vcf}"
    """
}

// BCFtools filter ~ extract all deletion and duplication records that pass the duphold false positive filter
process extractFpFilterPassingSvCalls_bcftools {
    tag "${tumor_normal_sample_id}"

    input:
    tuple val(tumor_normal_sample_id), path(consensus_somatic_sv_fpmarked_vcf) from consensus_sv_vcf_forFpFiltering

    output:
    tuple val(tumor_normal_sample_id), path(consensus_somatic_sv_fpfiltered_vcf) into consensus_sv_vcf_forAnnotation

    when:
    params.manta == "on" && params.svaba == "on" && params.delly == "on"

    script:
    consensus_somatic_sv_fpfiltered_vcf = "${tumor_normal_sample_id}.consensus.somatic.sv.fpfiltered.vcf"
    """
    bcftools filter \
    --output-type v \
    --exclude 'INFO/SVTYPE="DEL" && FORMAT/DHFFC>0.7' \
    "${consensus_somatic_sv_fpmarked_vcf}" \
    | \
    bcftools filter \
    --output-type v \
    --exclude 'INFO/SVTYPE="DUP" && FORMAT/DHBFC<1.3' \
    --output "${consensus_somatic_sv_fpfiltered_vcf}"
    """
}

// Depending on whether the reference files used for AnnotSV annotation was pre-downloaded, set the input
// channel for the AnnotSV SV annotation process
if( params.annotsv_ref_cached == "yes" ) {
    annotsv_ref_dir_forSvAnnotation = annotsv_ref_dir_pre_downloaded_forSvAnnotation
}
else {
    annotsv_ref_dir_forSvAnnotation = annotsv_ref_dir_from_process_forSvAnnotation
}

// AnnotSV ~ annotate consensus SV calls with multiple resources
process annotateConsensusSvCalls_annotsv {
    publishDir "${params.output_dir}/somatic/consensus/${tumor_normal_sample_id}", mode: 'copy', pattern: '*.{sv.annotated.genesplit.bed,bedpe}'
    tag  "${tumor_normal_sample_id}"

    input:
    tuple val(tumor_normal_sample_id), path(consensus_somatic_sv_fpfiltered_vcf), path(final_manta_somatic_sv_read_support), path(final_svaba_somatic_sv_read_support), path(final_delly_somatic_sv_read_support), path(annotsv_ref_dir_bundle) from consensus_sv_vcf_forAnnotation.join(manta_sv_read_support_forAnnotation).join(svaba_sv_read_support_forAnnotation).join(delly_sv_read_support_forAnnotation).combine(annotsv_ref_dir_forSvAnnotation)

    output:
    path gene_split_annotated_consensus_sv_bed
    path hq_consensus_sv_bedpe

    when:
    params.manta == "on" && params.svaba == "on" && params.delly == "on"

    script:
    gene_split_annotated_consensus_sv_bed = "${tumor_normal_sample_id}.hq.consensus.somatic.sv.annotated.genesplit.bed"
    collapsed_annotated_consensus_sv_bed = "${tumor_normal_sample_id}.hq.consensus.somatic.sv.annotated.collapsed.bed"
    hq_consensus_sv_bedpe = "${tumor_normal_sample_id}.hq.consensus.somatic.sv.annotated.bedpe"
    """
    grep -v 'SVTYPE=BND' "${consensus_somatic_sv_fpfiltered_vcf}" > "${tumor_normal_sample_id}.consensus.somatic.sv.nonbreakend.vcf"

	\$ANNOTSV/bin/AnnotSV \
	-annotationsDir "${annotsv_ref_dir_bundle}" \
	-annotationMode split \
	-genomeBuild GRCh38 \
	-hpo HP:0006775 \
	-outputDir . \
	-outputFile "${tumor_normal_sample_id}.consensus.somatic.sv.nonbreakend.annotated.genesplit" \
	-SVinputFile "${tumor_normal_sample_id}.consensus.somatic.sv.nonbreakend.vcf" \
	-SVminSize 1 \
	-includeCI 0 \
	-tx ENSEMBL

    paste \
	<(cut -f 2 "${tumor_normal_sample_id}.consensus.somatic.sv.nonbreakend.annotated.genesplit.tsv" | awk 'BEGIN {OFS="\t"} {print "chr"\$1}') \
	<(cut -f 3-4 "${tumor_normal_sample_id}.consensus.somatic.sv.nonbreakend.annotated.genesplit.tsv" | awk 'BEGIN {OFS="\t"} {print \$1-1,\$2}') \
	<(cut -f 20 "${tumor_normal_sample_id}.consensus.somatic.sv.nonbreakend.annotated.genesplit.tsv") \
	<(cut -f 5-6 "${tumor_normal_sample_id}.consensus.somatic.sv.nonbreakend.annotated.genesplit.tsv") \
	<(cut -f 1 "${tumor_normal_sample_id}.consensus.somatic.sv.nonbreakend.annotated.genesplit.tsv") \
	<(cut -f 21,23-36,64-65 "${tumor_normal_sample_id}.consensus.somatic.sv.nonbreakend.annotated.genesplit.tsv" | sed 's|\t\t|\t.\t|g' | sed 's|\t\t|\t.\t|g' | sed 's|\t\$|\t.|') \
	<(cut -f 44-45,48-49,52-63 "${tumor_normal_sample_id}.consensus.somatic.sv.nonbreakend.annotated.genesplit.tsv" | sed 's|^\t|.\t|' | sed 's|\t\t|\t.\t|g' | sed 's|\t\t|\t.\t|g') \
	| \
	sed 's|\t\$|\t.|' \
	| \
	sed 's|chrSV_chrom\t-1|SV_chrom\tSV_start|' \
	| \
	sort -k1,1V -k2,2n > "${tumor_normal_sample_id}.hq.consensus.somatic.sv.nonbreakend.annotated.genesplit.bed"

	\$ANNOTSV/bin/AnnotSV \
	-annotationsDir "${annotsv_ref_dir_bundle}" \
	-annotationMode full \
	-genomeBuild GRCh38 \
	-hpo HP:0006775 \
	-outputDir . \
	-outputFile "${tumor_normal_sample_id}.consensus.somatic.sv.nonbreakend.annotated.collapsed" \
	-SVinputFile "${tumor_normal_sample_id}.consensus.somatic.sv.nonbreakend.vcf" \
	-SVminSize 1 \
	-includeCI 0 \
	-tx ENSEMBL

	paste \
	<(cut -f 2 "${tumor_normal_sample_id}.consensus.somatic.sv.nonbreakend.annotated.collapsed.tsv" | awk 'BEGIN {OFS="\t"} {print "chr"\$1}') \
	<(cut -f 3-4 "${tumor_normal_sample_id}.consensus.somatic.sv.nonbreakend.annotated.collapsed.tsv") \
	<(cut -f 20 "${tumor_normal_sample_id}.consensus.somatic.sv.nonbreakend.annotated.collapsed.tsv") \
	<(cut -f 5-6 "${tumor_normal_sample_id}.consensus.somatic.sv.nonbreakend.annotated.collapsed.tsv") \
	<(cut -f 1,106,108 "${tumor_normal_sample_id}.consensus.somatic.sv.nonbreakend.annotated.collapsed.tsv" | sed 's|\t\t|\t.\t|g') \
	<(cut -f 8,13,15-17,21-22,66-79 "${tumor_normal_sample_id}.consensus.somatic.sv.nonbreakend.annotated.collapsed.tsv" | sed 's|\t\t|\t.\t|g' | sed 's|\t\t|\t.\t|g') \
	| \
	sed 's|\t\$|\t.|' \
	| \
	sed 's|chrSV_chrom|SV_chrom|' \
	| \
	sort -k1,1V -k2,2n > "${tumor_normal_sample_id}.hq.consensus.somatic.sv.nonbreakend.annotated.collapsed.bed"


	grep -E '^#|SVTYPE=BND' "${consensus_somatic_sv_fpfiltered_vcf}" > "${tumor_normal_sample_id}.consensus.somatic.sv.breakend.vcf"

	\$ANNOTSV/bin/AnnotSV \
	-annotationsDir "${annotsv_ref_dir_bundle}" \
	-annotationMode split \
	-genomeBuild GRCh38 \
	-hpo HP:0006775 \
	-outputDir . \
	-outputFile "${tumor_normal_sample_id}.consensus.somatic.sv.breakend.annotated.genesplit" \
	-SVinputFile "${tumor_normal_sample_id}.consensus.somatic.sv.breakend.vcf" \
	-SVminSize 1 \
	-includeCI 0 \
	-tx ENSEMBL

	paste \
    <(cut -f 2 "${tumor_normal_sample_id}.consensus.somatic.sv.breakend.annotated.genesplit.tsv" | awk 'BEGIN {OFS="\t"} {print "chr"\$1}') \
    <(cut -f 3-4 "${tumor_normal_sample_id}.consensus.somatic.sv.breakend.annotated.genesplit.tsv" | awk 'BEGIN {OFS="\t"} {print \$1-1,\$2}') \
    <(cut -f 20 "${tumor_normal_sample_id}.consensus.somatic.sv.breakend.annotated.genesplit.tsv") \
    <(cut -f 5-6 "${tumor_normal_sample_id}.consensus.somatic.sv.breakend.annotated.genesplit.tsv") \
    <(cut -f 1 "${tumor_normal_sample_id}.consensus.somatic.sv.breakend.annotated.genesplit.tsv") \
    <(cut -f 21,23-36,64-65 "${tumor_normal_sample_id}.consensus.somatic.sv.breakend.annotated.genesplit.tsv" | sed 's|\t\t|\t.\t|g' | sed 's|\t\t|\t.\t|g' | sed 's|\t\$|\t.|') \
    <(cut -f 44-45,48-49,52-63 "${tumor_normal_sample_id}.consensus.somatic.sv.breakend.annotated.genesplit.tsv" | sed 's|^\t|.\t|' | sed 's|\t\t|\t.\t|g' | sed 's|\t\t|\t.\t|g') \
    | \
    sed 's|\t\$|\t.|' \
    | \
    sed 's|chrSV_chrom\t-1|SV_chrom\tSV_start|' \
    | \
    sort -k1,1V -k2,2n > "${tumor_normal_sample_id}.hq.consensus.somatic.sv.breakend.annotated.genesplit.bed"

    \$ANNOTSV/bin/AnnotSV \
    -annotationsDir "${annotsv_ref_dir_bundle}" \
    -annotationMode full \
    -genomeBuild GRCh38 \
    -hpo HP:0006775 \
    -outputDir . \
    -outputFile "${tumor_normal_sample_id}.consensus.somatic.sv.breakend.annotated.collapsed" \
    -SVinputFile "${tumor_normal_sample_id}.consensus.somatic.sv.breakend.vcf" \
    -SVminSize 1 \
	-includeCI 0 \
    -tx ENSEMBL

    paste \
    <(cut -f 2 "${tumor_normal_sample_id}.consensus.somatic.sv.breakend.annotated.collapsed.tsv" | awk 'BEGIN {OFS="\t"} {print "chr"\$1}') \
    <(cut -f 3-4 "${tumor_normal_sample_id}.consensus.somatic.sv.breakend.annotated.collapsed.tsv") \
    <(cut -f 20 "${tumor_normal_sample_id}.consensus.somatic.sv.breakend.annotated.collapsed.tsv") \
    <(cut -f 5-6 "${tumor_normal_sample_id}.consensus.somatic.sv.breakend.annotated.collapsed.tsv") \
    <(cut -f 1,106,108 "${tumor_normal_sample_id}.consensus.somatic.sv.breakend.annotated.collapsed.tsv" | sed 's|\t\t|\t.\t|g') \
    <(cut -f 8,13,15-17,21-22,66-79 "${tumor_normal_sample_id}.consensus.somatic.sv.breakend.annotated.collapsed.tsv" | sed 's|\t\t|\t.\t|g' | sed 's|\t\t|\t.\t|g') \
    | \
    sed 's|\t\$|\t.|' \
    | \
    sed 's|chrSV_chrom|SV_chrom|' \
    | \
    sort -k1,1V -k2,2n > "${tumor_normal_sample_id}.hq.consensus.somatic.sv.breakend.annotated.collapsed.bed"

    cat "${tumor_normal_sample_id}.hq.consensus.somatic.sv.nonbreakend.annotated.genesplit.bed" \
    <(grep -v 'SV_chrom' "${tumor_normal_sample_id}.hq.consensus.somatic.sv.breakend.annotated.genesplit.bed") \
    | \
    sort -k1,1V -k2,2n > "${gene_split_annotated_consensus_sv_bed}"

    cat "${tumor_normal_sample_id}.hq.consensus.somatic.sv.nonbreakend.annotated.collapsed.bed" \
    <(grep -v 'SV_chrom' "${tumor_normal_sample_id}.hq.consensus.somatic.sv.breakend.annotated.collapsed.bed") \
    | \
    sort -k1,1V -k2,2n > "${collapsed_annotated_consensus_sv_bed}"


    # Transform collapsed annotation BED to high quality BEDPE
    high_quality_bedpe_transformer.py \
    <(grep -v 'SV_chrom' "${collapsed_annotated_consensus_sv_bed}") \
    "${final_manta_somatic_sv_read_support}" \
    "${final_svaba_somatic_sv_read_support}" \
    "${final_delly_somatic_sv_read_support}" \
    "${hq_consensus_sv_bedpe}"
    """
}

// END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \\


// ~~~~~~~~~~ CONSENSUS METADATA ~~~~~~~~~~~~ \\
// START

// Merge various metadata output, simple concatenation of per tool output
process mergeMetadataOutput {
	publishDir "${params.output_dir}/somatic/consensus/${tumor_normal_sample_id}", mode: 'copy', pattern: '*.{txt}'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(sample_sex), path(conpair_concordance_file), path(conpair_contamination_file), path(mutect_contamination_file), path(control_freec_run_info), path(sclust_cnv_summary_file) from allelecount_output_forConsensusMetadata.join(conpair_output_forConsensusMetadata).join(mutect_output_forConsensusMetadata).join(control_freec_output_forConsensusMetadata).join(sclust_output_forConsensusMetadata)

	output:
	path consensus_metadata_file

	when:
	params.conpair == "on" && params.mutect == "on" && params.sclust == "on"

	script:
	consensus_metadata_file = "${tumor_normal_sample_id}.consensus.somatic.metadata.txt"
	"""
	touch "${consensus_metadata_file}"

	echo "### alleleCount ###" >> "${consensus_metadata_file}"
	echo "" >> "${consensus_metadata_file}"
	echo "sample_sex" >> "${consensus_metadata_file}"
	cat "${sample_sex}" >> "${consensus_metadata_file}"
	echo "" >> "${consensus_metadata_file}"

	echo "### Conpair ###" >> "${consensus_metadata_file}"
	echo "" >> "${consensus_metadata_file}"
	cat "${conpair_concordance_file}" >> "${consensus_metadata_file}"
	echo "" >> "${consensus_metadata_file}"
	echo "cross-sample"
	cat "${conpair_contamination_file}" >> "${consensus_metadata_file}"
	echo "" >> "${consensus_metadata_file}"

	echo "### Mutect2 ###" >> "${consensus_metadata_file}"
	echo "" >> "${consensus_metadata_file}"
	echo "cross-sample"
	cut -f 2,3 "${mutect_contamination_file}" >> "${consensus_metadata_file}"
	echo "" >> "${consensus_metadata_file}"

	echo "### Control-FREEC ###" >> "${consensus_metadata_file}"
	echo "" >> "${consensus_metadata_file}"
	grep 'Sample_Purity' "${control_freec_run_info}" >> "${consensus_metadata_file}"
	grep 'Output_Ploidy' "${control_freec_run_info}" >> "${consensus_metadata_file}"
	echo "" >> "${consensus_metadata_file}"

	echo "### Sclust ###" >> "${consensus_metadata_file}"
	echo "" >> "${consensus_metadata_file}"
	cut -f 2,3,5 "${sclust_cnv_summary_file}" >> "${consensus_metadata_file}"
	echo "" >> "${consensus_metadata_file}"
	"""
}

// END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \\


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

// VEP ~ annotate the final somatic SNV and InDel VCFs using databases including Ensembl, GENCODE, RefSeq, PolyPhen, SIFT, dbSNP, COSMIC, etc.
process annotateSnvAndIndelVcf_vep {
	publishDir "${params.output_dir}/somatic/consensus/${tumor_normal_sample_id}", mode: 'copy', pattern: '*.{vcf.gz,html}'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(hq_snv_consensus_vcf), path(hq_snv_consensus_vcf_index), path(hq_indel_consensus_vcf), path(hq_indel_consensus_vcf_index), path(cached_ref_dir_vep), path(reference_genome_fasta_forAnnotation), path(reference_genome_fasta_index_forAnnotation), path(reference_genome_fasta_dict_forAnnotation) from high_quality_consensus_snv_forAnnotation.join(high_quality_consensus_indel_forAnnotation).combine(vep_ref_dir).combine(reference_genome_bundle_forAnnotation)

	output:
	path("*.annotated.vcf.gz")
	path("*.vep.summary.html")

	when:
	params.varscan == "on" && params.mutect == "on" && params.strelka == "on" && params.svaba == "on"

	script:
	"""
	for vcf in `ls *.vcf.gz`;
		do
			output_vcf=\$(echo \${vcf} | sed 's|.vcf.gz|.annotated.vcf.gz|')
			output_stats=\$(echo \${vcf} | sed 's|.vcf.gz|.vep.summary.html|')

			vep \
			--offline \
			--cache \
			--dir "${cached_ref_dir_vep}" \
			--assembly GRCh38 \
			--fasta "${reference_genome_fasta_forAnnotation}" \
			--input_file \${vcf} \
			--format vcf \
			--sift b \
			--polyphen b \
			--biotype \
			--symbol \
			--nearest symbol \
			--check_existing \
			--max_af \
			--gencode_basic \
			--pick \
			--stats_file \${output_stats} \
			--output_file \${output_vcf} \
			--compress_output bgzip \
			--vcf
		done
	"""
}

// END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \\
