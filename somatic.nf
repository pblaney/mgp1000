// Myeloma Genome Project 1000
// Comprehensive pipeline for analysis of matched T/N Multiple Myeloma WGS data
// https://github.com/pblaney/mgp1000

// This portion of the pipeline is used for somatic variant analysis of matched tumor/normal WGS samples.
// It is designed to be run with BAMs that were genereated via the Preprocessing module of this pipeline.

import java.text.SimpleDateFormat;
def workflowTimestamp = "${workflow.start.format('MM-dd-yyyy HH:mm')}"

def helpMessage() {
	log.info"""
	                             .------------------------.
	                            |    .-..-. .--. .---.     |
	                            |    : `' :: .--': .; :    |
	                            |    : .. :: : _ :  _.'    |
	                            |    : :; :: :; :: :       |
	                            |    :_;:_;`.__.':_;       |
	                            |   ,-. .--.  .--.  .--.   |
	                            | .'  :: ,. :: ,. :: ,. :  |
	                            |   : :: :: :: :: :: :: :  |
	                            |   : :: :; :: :; :: :; :  |
	                            |   :_;`.__.'`.__.'`.__.'  |
	                             .________________________.

	                            ░█▀▀░█▀█░█▄█░█▀█░▀█▀░▀█▀░█▀▀
                                ░▀▀█░█░█░█░█░█▀█░░█░░░█░░█░░
                                ░▀▀▀░▀▀▀░▀░▀░▀░▀░░▀░░▀▀▀░▀▀▀

	Usage:
	  nextflow run somatic.nf --run_id STR --sample_sheet FILE -profile somatic [-bg] [-resume]
	  [--input_dir PATH] [--output_dir PATH] [--email STR] [--mutect_ref_vcf_concatenated STR]
	  [--battenberg_ref_cached STR] [--annotsv_ref_cached STR] [--vep_ref_cached STR] [--help]

	Mandatory Arguments:
	  --run_id                       STR  Unique identifier for pipeline run
	  --sample_sheet                FILE  CSV file containing the list of samples where the
	                                      first column designates the file name of the normal
	                                      sample, the second column for the file name of the
	                                      matched tumor sample
	  -profile                       STR  Configuration profile to use, must use somatic

	Main Options:
	  -bg                           FLAG  Runs the pipeline processes in the background, this
	                                      option should be included if deploying pipeline with
	                                      real data set so processes will not be cut if user
	                                      disconnects from deployment environment
	  -resume                       FLAG  Successfully completed tasks are cached so that if
	                                      the pipeline stops prematurely the previously
	                                      completed tasks are skipped while maintaining their
	                                      output
	  --input_dir                   PATH  Directory that holds BAMs and associated index files,
	                                      this should be given as an absolute path
	                                      [Default: input/preprocessedBams/]
	  --output_dir                  PATH  Directory that will hold all output files this should
	                                      be given as an absolute path
	                                      [Default: output/]
	  --email                        STR  Email address to send workflow completion/stoppage
	                                      notification
	  --mutect_ref_vcf_concatenated  STR  Indicates whether or not the gnomAD allele frequency
	                                      reference VCF used for MuTect2 processes has been
	                                      concatenated, this will be done in a process of the
	                                      pipeline if it has not, this does not need to be done
	                                      for every separate run after the first
	                                      [Default: no | Available: yes, no]
	  --battenberg_ref_cached        STR  Indicates whether or not the reference files used for
	                                      Battenberg have been downloaded/cached locally, this
	                                      will be done in a process of the pipeline if it has
	                                      not, this does not need to be done for every separate
	                                      run after the first
	                                      [Default: no | Available: yes, no]
	  --cpus                         INT  Globally set the number of cpus to be allocated
	  --memory                       STR  Globally set the amount of memory to be allocated,
	                                      written as '##.GB' or '##.MB'
	  --queue_size                   INT  Set max number of tasks the pipeline will launch
	                                      [Default: 100]
	  --executor                     STR  Set the job executor for the run
	                                      [Default: slurm | Available: local, slurm, lsf]
	  --help                        FLAG  Prints this message

	Toolbox Options:
	  --battenberg                   STR  Indicates whether or not to use this tool
	                                      [Default: on | Available: off, on]
	  --battenberg_min_depth         STR  Manually set the minimum read depth in the normal
	                                      sample for SNP filtering in BAF calculations,
	                                      default is for 30x coverage
	                                      [Default: 10]
	  --battenberg_preset_rho_psi    STR  Wish to manually set the rho/psi for this run?
	                                      If TRUE, must set both rho and psi.
	                                      [Default: FALSE | Available: FALSE, TRUE]
	  --battenberg_preset_rho        INT  Manually set the value of rho (purity)
	                                      [Default: NA]
	  --battenberg_preset_psi        INT  Manually set the value of psi (ploidy)
	                                      [Default: NA]
	  --facets                       STR  Indicates whether or not to use this tool
	                                      [Default: on | Available: off, on]
	  --facets_min_depth             INT  Manually set the minimum read depth in the normal
	                                      sample for SNP filtering in BAF calculations,
	                                      default is for 30x coverage
	                                      [Default: 20]
	  --manta                        STR  Indicates whether or not to use this tool
	                                      [Default: on | Available: off, on]
	  --svaba                        STR  Indicates whether or not to use this tool
	                                      [Default: on | Available: off, on]
	  --delly                        STR  Indicates whether or not to use this tool
	                                      [Default: on | Available: off, on]
	  --delly_strict                 STR  Enforce stricter thresholds for calling SVs with DELLY
	                                      to overcome libraries with extraordinary number of
	                                      interchromosomal reads
	                                      [Default: off | Available: off, on]
	  --igcaller                     STR  Indicates whether or not to use this tool
	                                      [Default: on | Available: off, on]
	  --varscan                      STR  Indicates whether or not to use this tool
	                                      [Default: on | Available: off, on]
	  --mutect                       STR  Indicates whether or not to use this tool
	                                      [Default: on | Available: off, on]
	  --strelka                      STR  Indicates whether or not to use this tool
	                                      [Default: on | Available: off, on]
	  --conpair                      STR  Indicates whether or not to use this tool
	                                      [Default: on | Available: off, on]
	  --conpair_min_cov              INT  Manually set the minimum coverage
	                                      [Default: 10]
	  --fragcounter                  STR  Indicates whether or not to use this tool
	                                      [Default: on | Available: off, on]
	  --telomerecat                  STR  Indicates whether or not to use this tool
	                                      [Default: off | Available: off, on]
	  --telomerehunter               STR  Indicates whether or not to use this tool
	                                      [Default: off | Available: off, on]

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
params.seq_protocol = "WGS"
params.mutect_ref_vcf_concatenated = "yes"
params.battenberg_ref_cached = "yes"
params.vep_ref_cached = "yes"
params.telomerecat = "off"
params.telomerehunter = "off"
params.conpair = "on"
params.varscan = "on"
params.mutect = "on"
params.strelka = "on"
params.caveman = "off"
params.fragcounter = "on"
params.battenberg = "on"
params.facets = "on"
params.manta = "on"
params.svaba = "on"
params.delly = "on"
params.delly_strict = "off"
params.igcaller = "on"
params.conpair_min_cov = 10
params.battenberg_min_depth = 10
params.battenberg_preset_rho_psi = "FALSE"
params.battenberg_preset_rho = "NA"
params.battenberg_preset_psi = "NA"
params.facets_min_depth = 20
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

// Print preemptive error message if Strelka is set while Manta is not
if( params.strelka == "on" && params.manta == "off" ) exit 1, "Strelka requires output from Manta to run so both must be turned on"

// Set channels for reference files
Channel
    .value( file('references/hg38/Homo_sapiens_assembly38.fasta') )
    .set{ ref_genome_fasta_file }

Channel
    .value( file('references/hg38/Homo_sapiens_assembly38.fasta.fai') )
    .set{ ref_genome_fasta_index_file }

Channel
    .value( file('references/hg38/Homo_sapiens_assembly38.dict') )
    .set{ ref_genome_fasta_dict_file }

Channel
    .value( file('references/hg38/wgs_calling_regions.hg38.bed') )
    .set{ gatk_wgs_bed }

Channel
    .value( file('references/hg38/wgs_calling_regions_blacklist.0based.hg38.bed') )
    .set{ wgs_bed_blacklist_0based }

Channel
	.fromList( ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6',
	            'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
	            'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18',
	            'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY',] )
	.into{ chromosome_list_forVarscanSamtoolsMpileup;
	       chromosome_list_forMutectCalling;
	       chromosome_list_forMutectPileup;
	       chromosome_list_forCavemanSplit }

Channel
    .value( file('references/hg38/sex_identification_loci.chrY.hg38.txt') )
    .set{ sex_identification_loci }

Channel
    .value( file('references/hg38/cytoband_autosome_sex_chroms.hg38.bed') )
    .set{ cytoband_bed }

Channel
    .value( file('references/hg38/1000g_pon.hg38.vcf.gz') )
    .set{ panel_of_normals_1000G }

Channel
    .value( file('references/hg38/1000g_pon.hg38.vcf.gz.tbi') )
    .set{ panel_of_normals_1000G_index }	

Channel
    .value( file('references/hg38/af-only-gnomad.chr1-9.hg38.vcf.gz') )
    .set{ gnomad_ref_vcf_chromosomes1_9 }

Channel
    .value( file('references/hg38/af-only-gnomad.chr1-9.hg38.vcf.gz.tbi') )
    .set{ gnomad_ref_vcf_chromosomes1_9_index }

Channel
    .value( file('references/hg38/af-only-gnomad.chr10-22.hg38.vcf.gz') )
    .set{ gnomad_ref_vcf_chromosomes10_22 }

Channel
    .value( file('references/hg38/af-only-gnomad.chr10-22.hg38.vcf.gz.tbi') )
    .set{ gnomad_ref_vcf_chromosomes10_22_index }

Channel
    .value( file('references/hg38/af-only-gnomad.chrXYM-alts.hg38.vcf.gz') )
    .set{ gnomad_ref_vcf_chromosomesXYM_alts }

Channel
    .value( file('references/hg38/af-only-gnomad.chrXYM-alts.hg38.vcf.gz.tbi') )
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
    .value( file('references/hg38/small_exac_common_3.hg38.vcf.gz') )
    .set{ exac_common_sites_ref_vcf }

Channel
    .value( file('references/hg38/small_exac_common_3.hg38.vcf.gz.tbi') )
    .set{ exac_common_sites_ref_vcf_index }

if( params.battenberg_ref_cached == "yes" ) {
	Channel
		.fromPath( 'references/hg38/battenberg_reference/', checkIfExists: true )
		.ifEmpty{ error "The run command issued has the '--battenberg_ref_cached' parameter set to 'yes', however the directory does not exist. Please set the '--battenberg_ref_cached' parameter to 'no' and resubmit the run command. For more information, check the README or issue the command 'nextflow run somatic.nf --help'"}
		.set{ battenberg_ref_dir_preDownloaded }
}

Channel
    .value( file('references/hg38') )
    .set{ gc_mappability_dir }

Channel
    .value( file('references/hg38/Hapmap_3.3.hg38.vcf.gz') )
    .set{ hapmap_ref_snps_vcf }

Channel
    .value( file('references/hg38/Hapmap_3.3.hg38.vcf.gz.tbi') )
    .set{ hapmap_ref_snps_vcf_index }

Channel
    .value( file('references/hg38/common_all_20180418.vcf.gz') )
    .set{ common_dbsnp_ref_vcf }

Channel
    .value( file('references/hg38/common_all_20180418.vcf.gz.tbi') )
    .set{ common_dbsnp_ref_vcf_index }

Channel
	.fromPath( ['references/hg38/Homo_sapiens_assembly38.fasta', 'references/hg38/Homo_sapiens_assembly38.fasta.fai',
	            'references/hg38/Homo_sapiens_assembly38.fasta.64.alt', 'references/hg38/Homo_sapiens_assembly38.fasta.64.amb',
	            'references/hg38/Homo_sapiens_assembly38.fasta.64.ann', 'references/hg38/Homo_sapiens_assembly38.fasta.64.bwt',
	            'references/hg38/Homo_sapiens_assembly38.fasta.64.pac', 'references/hg38/Homo_sapiens_assembly38.fasta.64.sa'] )
	.set{ bwa_ref_genome_files }

Channel
    .value( file('references/hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz') )
    .set{ dbsnp_known_indel_ref_vcf }

Channel
    .value( file('references/hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi') )
    .set{ dbsnp_known_indel_ref_vcf_index }

Channel
    .value( file('references/hg38/simple_and_centromeric_repeats.hg38.bed') )
    .set{ simple_and_centromeric_repeats_bed }

Channel
	.value( file('references/hg38/wgs_calling_regions_blacklist.1based.hg38.bed') )
	.set{ wgs_bed_blacklist_0based }

Channel
	.value( file('references/hg38/exome_blacklist.hg38.bed') )
	.set{ wes_bed_blacklist_0based }	

Channel
	.value( file('references/hg38/unmatchedNormal.bed.gz') )
	.set{ unmatched_normal_bed }

Channel
	.value( file('references/hg38/unmatchedNormal.bed.gz.tbi') )
	.set{ unmatched_normal_bed_index }

Channel
	.value( file('references/hg38/centromeric_repeats.hg38.bed.gz') )
	.set{ centromeric_repeats_bed }

Channel
	.value( file('references/hg38/centromeric_repeats.hg38.bed.gz.tbi') )
	.set{ centromeric_repeats_bed_index }

Channel
	.value( file('references/hg38/simple_repeats.hg38.bed.gz') )
	.set{ simple_repeats_bed }

Channel
	.value( file('references/hg38/simple_repeats.hg38.bed.gz.tbi') )
	.set{ simple_repeats_bed_index }

Channel
	.value( file('references/hg38/dbsnp138.hg38.bed.gz') )
	.set{ dbsnp_bed }

Channel
	.value( file('references/hg38/dbsnp138.hg38.bed.gz.tbi') )
	.set{ dbsnp_bed_index }

Channel
    .value( file('references/hg38/Homo_sapiens.GRCh38.108.gtf.gz') )
    .set{ ensembl_gtf_file }


// #################################################### \\
// ~~~~~~~~~~~~~~~~ PIPELINE PROCESSES ~~~~~~~~~~~~~~~~ \\

log.info ''
log.info '################################################'
log.info ''
log.info "           .------------------------.           "
log.info "          |    .-..-. .--. .---.     |          "
log.info "          |    : `' :: .--': .; :    |          "
log.info "          |    : .. :: : _ :  _.'    |          "
log.info "          |    : :; :: :; :: :       |          "
log.info "          |    :_;:_;`.__.':_;       |          "
log.info "          |   ,-. .--.  .--.  .--.   |          "
log.info "          | .'  :: ,. :: ,. :: ,. :  |          "
log.info "          |   : :: :: :: :: :: :: :  |          "
log.info "          |   : :: :; :: :; :: :; :  |          "
log.info "          |   :_;`.__.'`.__.'`.__.'  |          "
log.info "           .________________________.           "
log.info ''
log.info "          ░█▀▀░█▀█░█▄█░█▀█░▀█▀░▀█▀░█▀▀           "
log.info "          ░▀▀█░█░█░█░█░█▀█░░█░░░█░░█░░           "
log.info "          ░▀▀▀░▀▀▀░▀░▀░▀░▀░░▀░░▀▀▀░▀▀▀           "
log.info ''
log.info "~~~ Launch Time ~~~"
log.info ''
log.info " ${workflowTimestamp}"
log.info ''
log.info "~~~ Input Directory ~~~"
log.info ''
log.info " ${params.input_dir}"
log.info ''
log.info "~~~ Output Directory ~~~"
log.info ''
log.info " ${params.output_dir}"
log.info ''
log.info "~~~ Run Report File ~~~"
log.info ''
log.info " nextflow_report.${params.run_id}.html"
log.info ''
log.info "~~~ Sequencing Protocol ~~~"
log.info ''
log.info " ${params.seq_protocol}"
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
	       tumor_normal_pair_forManta;
	       tumor_normal_pair_forSvaba;
	       tumor_normal_pair_forDelly;
	       tumor_normal_pair_forIgCaller }

// alleleCount ~ determine the sex of each sample to use in downstream analyses
process identifySampleSex_allelecount {
	publishDir "${params.output_dir}/somatic/sexOfSamples", mode: 'copy', pattern: '*.{txt}'
	tag "${tumor_normal_sample_id}"

	input:
	tuple path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index) from tumor_normal_pair_forAlleleCount
	path ref_genome_fasta_index from ref_genome_fasta_index_file
	path sex_loci from sex_identification_loci

	output:
	tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index) into bams_forVarscanBamReadcount
	tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index) into tumor_bams_forFragCounter
	tuple val(tumor_normal_sample_id), path(normal_bam), path(normal_bam_index) into normal_bams_forFragCounter
	tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index) into bams_forFragCounterPileup
	tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(sample_sex) into bams_and_sex_of_sample_forBattenberg
	tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index) into bams_forFacetsPileup
	tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index) into bams_forCaveman
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
	--loci-file "${sex_loci}" \
	--hts-file "${normal_bam}" \
	--ref-file "${ref_genome_fasta_index}" \
	--output-file "${sex_loci_allele_counts}"

	sample_sex_determinator.sh "${sex_loci_allele_counts}" > "${sample_sex}"
	"""
}


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
  	battenberg_references = "battenberg_reference/"
  	"""
  	mkdir -p ${battenberg_references}
  	cd ${battenberg_references}/
  	mkdir -p temp/
  	mkdir -p GC_correction_hg38/
  	mkdir -p RT_correction_hg38/
  	mkdir -p shapeit2/
  	mkdir -p imputation/
  	mkdir -p 1000G_loci_hg38/
  	mkdir -p probloci/
  	mkdir -p beagle5/

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
  	tuple val(tumor_normal_sample_id), path(battenberg_cnv_profile) into tumor_cnv_profile_forCaveman

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
  	cp /opt/downloads/beagle.08Feb22.fa4.jar battenberg_reference/beagle5/

  	sex=\$(cut -d ' ' -f 1 "${sample_sex}")

  	battenberg_executor.sh \
  	"${tumor_id}" \
  	"${normal_id}" \
  	"${tumor_bam}" \
  	"${normal_bam}" \
  	"\${sex}" \
  	"${output_dir}" \
  	"${task.cpus}" \
  	"${params.battenberg_min_depth}" \
  	${params.battenberg_preset_rho_psi} \
    ${params.battenberg_preset_rho} \
    ${params.battenberg_preset_psi}

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


// ~~~~~~~~~~~~~~~~~ FACETS ~~~~~~~~~~~~~~~~ \\
// START

// FACETS ~ generate SNP read count pileups for CNV calling
process snpPileup_facets {
    publishDir "${params.output_dir}/somatic/facets", mode: 'copy'
    tag "${tumor_normal_sample_id}"

    input:
    tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index) from bams_forFacetsPileup
    path common_dbsnp_vcf from common_dbsnp_ref_vcf
	path common_dbsnp_vcf_index from common_dbsnp_ref_vcf_index

    output:
    tuple val(tumor_normal_sample_id), path(facets_snp_pileup) into snp_pileup_forFacets

    when:
    params.facets == "on"

    script:
    tumor_id = "${tumor_bam.baseName}".replaceFirst(/\..*$/, "")
    normal_id = "${normal_bam.baseName}".replaceFirst(/\..*$/, "")
    tumor_normal_sample_id = "${tumor_id}_vs_${normal_id}"
    facets_snp_pileup = "${tumor_normal_sample_id}.facets.snp_pileup.csv.gz"
    """
    snp-pileup \
    --gzip \
    --min-map-quality 15 \
    --min-base-quality 20 \
    --pseudo-snps 100 \
    --min-read-counts ${params.facets_min_depth} \
    "${common_dbsnp_vcf}" \
    "${facets_snp_pileup}" \
    "${normal_bam}" \
    "${tumor_bam}"
    """
}

// FACETS ~ fraction and copy number estimate from tumor/normal sequencing
process cnvCalling_facets {
    publishDir "${params.output_dir}/somatic/facets", mode: 'copy'
    tag "${tumor_normal_sample_id}"

    input:
    tuple val(tumor_normal_sample_id), path(facets_snp_pileup) from snp_pileup_forFacets

    output:
    tuple val(tumor_normal_sample_id), path(facets_cnv_profile) into facets_cnv_profile_forConsensusPrep
    path facets_purity_ploidy
    path facets_run_log
    path facets_cnv_pdf
    path facets_spider_qc_pdf

    when:
    params.facets == "on"

    script:
    facets_run_log = "${tumor_normal_sample_id}.facets.log.txt"
    facets_purity_ploidy = "${tumor_normal_sample_id}.facets.purity.ploidy.txt"
    facets_cnv_profile = "${tumor_normal_sample_id}.facets.cnv.txt"
    facets_cnv_pdf = "${tumor_normal_sample_id}.facets.cnv.pdf"
    facets_spider_qc_pdf = "${tumor_normal_sample_id}.facets.spider.pdf"
    """
    Rscript --vanilla ${workflow.projectDir}/bin/run_iarc_facets.R \
    "${facets_snp_pileup}" \
    hg38 \
    1000 \
    35 \
    150 \
    300 \
    ${params.facets_min_depth}

    mv "${tumor_normal_sample_id}.facets.snp_pileup.R_sessionInfo.txt" "${facets_run_log}"
    mv "${tumor_normal_sample_id}.facets.snp_pileup.def_cval"*"_stats.txt" "${facets_purity_ploidy}"
    mv "${tumor_normal_sample_id}.facets.snp_pileup.def_cval"*"_CNV.txt" "${facets_cnv_profile}"
    mv "${tumor_normal_sample_id}.facets.snp_pileup.def_cval"*"_CNV.pdf" "${facets_cnv_pdf}"
    mv "${tumor_normal_sample_id}.facets.snp_pileup.def_cval"*"_CNV_spider.pdf" "${facets_spider_qc_pdf}"
    """
}

// FACETS ~ fraction and copy number estimate from tumor/normal sequencing
process consensusCnvPrep_facets {
    tag "${tumor_normal_sample_id}"

    input:
    tuple val(tumor_normal_sample_id), path(facets_cnv_profile) from facets_cnv_profile_forConsensusPrep

    output:
    tuple val(tumor_normal_sample_id), path(facets_somatic_cnv_bed), path(facets_somatic_alleles_bed) into final_facets_cnv_profile_forConsensus

    when:
    params.facets == "on"

    script:
    facets_somatic_cnv_bed = "${tumor_normal_sample_id}.facets.somatic.cnv.bed"
    facets_somatic_alleles_bed = "${tumor_normal_sample_id}.facets.somatic.alleles.bed"
    """
    facets_cnv_profile_postprocesser.sh \
    "${facets_cnv_profile}" \
    "${tumor_normal_sample_id}" \
    "${facets_somatic_cnv_bed}" \
    "${facets_somatic_alleles_bed}"
    """
}

// END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \\


// ~~~~~~~~~~~~~~~~~ Manta ~~~~~~~~~~~~~~~~~ \\
// START

// Manta ~ call structural variants and indels from mapped paired-end sequencing reads of matched tumor/normal sample pairs
process svAndIndelCalling_manta {
	publishDir "${params.output_dir}/somatic/manta", mode: 'copy', pattern: '*.{vcf.gz,tbi}'
	tag "${tumor_normal_sample_id}"

	input:
	tuple path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index) from tumor_normal_pair_forManta
	path ref_genome_fasta from ref_genome_fasta_file
	path ref_genome_fasta_index from ref_genome_fasta_index_file
	path ref_genome_fasta_dict from ref_genome_fasta_dict_file
	path wgs_bed from gatk_wgs_bed

	output:
	tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(candidate_indel_vcf), path(candidate_indel_vcf_index) into bams_and_candidate_indel_vcf_forStrelka
	tuple val(tumor_normal_sample_id), val(tumor_id), val(normal_id), path(manta_somatic_sv_vcf), path(manta_somatic_sv_vcf_index) into manta_sv_vcf_forPostprocessing
	tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index) into manta_tumor_bam_forDuphold
	tuple path(unfiltered_sv_vcf), path(unfiltered_sv_vcf_index)
	tuple path(germline_sv_vcf), path(germline_sv_vcf_index)
	tuple val(tumor_normal_sample_id), path(germline_sv_vcf), path(germline_sv_vcf_index) into germline_indel_vcf_forCaveman

	when:
	params.manta == "on"

	script:
	tumor_id = "${tumor_bam.baseName}".replaceFirst(/\..*$/, "")
	normal_id = "${normal_bam.baseName}".replaceFirst(/\..*$/, "")
	tumor_normal_sample_id = "${tumor_id}_vs_${normal_id}"
	zipped_wgs_bed_forManta = "${wgs_bed}.gz"
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
	bgzip < "${wgs_bed}" > "${zipped_wgs_bed_forManta}"
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
	--referenceFasta "${ref_genome_fasta}" \
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
	"${ref_genome_fasta}" \
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
    tuple val(tumor_normal_sample_id), path(final_manta_somatic_sv_vcf) into manta_sv_vcf_forDuphold
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
    --output "${final_manta_somatic_sv_vcf}" \
    "${manta_somatic_sv_vcf}" 

    bcftools query \
    --format '%ID\t[%PR{1}]\t[%SR{1}]\n' \
    --output "${final_manta_somatic_sv_read_support}" \
    "${final_manta_somatic_sv_vcf}"
    """
}

// duphold ~ efficiently annotate SV calls with sequence depth information to reduce false positive deletion and duplication calls
process falsePostiveSvFilteringManta_duphold {
    publishDir "${params.output_dir}/somatic/manta", mode: 'copy', pattern: '*.{vcf.gz}'
    tag "${tumor_normal_sample_id}"

    input:
    tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index), path(final_manta_somatic_sv_vcf) from manta_tumor_bam_forDuphold.join(manta_sv_vcf_forDuphold)
    path ref_genome_fasta from ref_genome_fasta_file
    path ref_genome_fasta_index from ref_genome_fasta_index_file
    path ref_genome_fasta_dict from ref_genome_fasta_dict_file

    output:
    tuple val(tumor_normal_sample_id), path(manta_filtered_final_sv_vcf) into manta_filtered_final_sv_vcf_forConsensus

    when:
    params.manta == "on"

    script:
    manta_filtered_final_sv_vcf = "${tumor_normal_sample_id}.manta.somatic.sv.final.vcf.gz"
    """
    duphold \
    --vcf "${final_manta_somatic_sv_vcf}" \
    --bam "${tumor_bam}" \
    --fasta "${ref_genome_fasta}" \
    --threads ${task.cpus} \
    --output "${tumor_normal_sample_id}.manta.somatic.sv.fpmarked.vcf"

    # Filter using recommended thresholds for DEL/DUP
    bcftools filter \
    --output-type u \
    --exclude 'INFO/SVTYPE="DEL" && FORMAT/DHFFC>0.7' \
    "${tumor_normal_sample_id}.manta.somatic.sv.fpmarked.vcf" \
    | \
    bcftools filter \
    --output-type z \
    --exclude 'INFO/SVTYPE="DUP" && FORMAT/DHBFC<1.3' \
    --output "${manta_filtered_final_sv_vcf}"
    """
}

// END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \\


// ~~~~~~~~~~~~~~~~ Strelka2 ~~~~~~~~~~~~~~~ \\
// START

// Strelka2 ~ call germline and somatic small variants from mapped sequencing reads
process snvAndIndelCalling_strelka {
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(candidate_indel_vcf), path(candidate_indel_vcf_index) from bams_and_candidate_indel_vcf_forStrelka
	path ref_genome_fasta from ref_genome_fasta_file
	path ref_genome_fasta_index from ref_genome_fasta_index_file
	path ref_genome_fasta_dict from ref_genome_fasta_dict_file
	path wgs_bed from gatk_wgs_bed

	output:
	tuple val(tumor_normal_sample_id), path(unfiltered_strelka_snv_vcf), path(unfiltered_strelka_snv_vcf_index), path(unfiltered_strelka_indel_vcf), path(unfiltered_strelka_indel_vcf_index) into unfiltered_snv_and_indel_vcfs_forStrelkaBcftools

	when:
	params.strelka == "on" && params.manta == "on"

	script:
	zipped_wgs_bed_forStrelka = "${wgs_bed}.gz"
	unfiltered_strelka_snv_vcf = "${tumor_normal_sample_id}.strelka.somatic.snv.unfiltered.vcf.gz"
	unfiltered_strelka_snv_vcf_index = "${unfiltered_strelka_snv_vcf}.tbi"
	unfiltered_strelka_indel_vcf = "${tumor_normal_sample_id}.strelka.somatic.indel.unfiltered.vcf.gz"
	unfiltered_strelka_indel_vcf_index = "${unfiltered_strelka_indel_vcf}.tbi"
	"""
	bgzip < "${wgs_bed}" > "${zipped_wgs_bed_forStrelka}"
	tabix "${zipped_wgs_bed_forStrelka}"

	python \${STRELKA_DIR}/bin/configureStrelkaSomaticWorkflow.py \
	--normalBam "${normal_bam}" \
	--tumorBam "${tumor_bam}" \
	--referenceFasta "${ref_genome_fasta}" \
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

// BCFtools norm ~ split multiallelic sites into multiple rows then left-align and normalize indels
process splitMultiallelicAndLeftNormalizeStrelkaVcf_bcftools {
	publishDir "${params.output_dir}/somatic/strelka", mode: 'copy', pattern: '*.{vcf.gz,tbi,txt}'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(unfiltered_strelka_snv_vcf), path(unfiltered_strelka_snv_vcf_index), path(unfiltered_strelka_indel_vcf), path(unfiltered_strelka_indel_vcf_index) from unfiltered_snv_and_indel_vcfs_forStrelkaBcftools
	path ref_genome_fasta from ref_genome_fasta_file
	path ref_genome_fasta_index from ref_genome_fasta_index_file
	path ref_genome_fasta_dict from ref_genome_fasta_dict_file

	output:
	tuple val(tumor_normal_sample_id), path(final_strelka_snv_vcf) into final_strelka_snv_vcf_forConsensus
	tuple val(tumor_normal_sample_id), path(final_strelka_indel_vcf) into final_strelka_indel_vcf_forConsensus
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
	--output-type u \
	- 2>"${strelka_indel_multiallelics_stats}" \
	| \
	bcftools norm \
	--threads ${task.cpus} \
	--fasta-ref "${ref_genome_fasta}" \
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

// Combine all needed reference FASTA files into one channel for use in SvABA process
//bwa_ref_genome_files.collect()
//	.set{ bwa_ref_genome_and_wgs_bed }

// SvABA ~ detecting structural variants using genome-wide local assembly
process svAndIndelCalling_svaba {
	publishDir "${params.output_dir}/somatic/svaba", mode: 'copy', pattern: '*.{somatic.sv.unprocessed.vcf.gz,somatic.sv.unprocessed.vcf.gz.tbi}'
	tag "${tumor_normal_sample_id}"

	input:
	tuple path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path("Homo_sapiens_assembly38.fasta"), path("Homo_sapiens_assembly38.fasta.fai"), path("Homo_sapiens_assembly38.fasta.64.alt"), path("Homo_sapiens_assembly38.fasta.64.amb"), path("Homo_sapiens_assembly38.fasta.64.ann"), path("Homo_sapiens_assembly38.fasta.64.bwt"), path("Homo_sapiens_assembly38.fasta.64.pac"), path("Homo_sapiens_assembly38.fasta.64.sa") from tumor_normal_pair_forSvaba.combine(bwa_ref_genome_files.collect())
	path ref_genome_fasta_dict from ref_genome_fasta_dict_file
	path wgs_blacklist_0based_bed from wgs_bed_blacklist_0based
	path dbsnp_known_indel_vcf from dbsnp_known_indel_ref_vcf
	path dbsnp_known_indel_vcf_index from dbsnp_known_indel_ref_vcf_index
	path simple_and_centromeric_repeats from simple_and_centromeric_repeats_bed

	output:
	tuple val(tumor_normal_sample_id), path(filtered_somatic_indel_vcf), path(filtered_somatic_indel_vcf_index) into filtered_indel_vcf_forSvabaBcftools
	tuple val(tumor_normal_sample_id), val(tumor_id), path(svaba_somatic_sv_vcf), path(svaba_somatic_sv_vcf_index), path(svaba_somatic_sv_unclassified_vcf), path(sample_renaming_file) into svaba_sv_vcf_forPostprocessing
	tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index) into svaba_tumor_bam_forDuphold
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
	--blacklist "${wgs_blacklist_0based_bed}" \
	--id-string "${tumor_normal_sample_id}" \
	--dbsnp-vcf "${dbsnp_known_indel_vcf}" \
	--simple-seq-database "${simple_and_centromeric_repeats_bed}" \
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
// read mapping quality, prepare VCF for consensus
process filterAndPostprocessSvabaVcf_bcftools {
    tag "${tumor_normal_sample_id}"

    input:
    tuple val(tumor_normal_sample_id), val(tumor_id), path(svaba_somatic_sv_vcf), path(svaba_somatic_sv_vcf_index), path(svaba_somatic_sv_unclassified_vcf), path(sample_renaming_file) from svaba_sv_vcf_forPostprocessing

    output:
    tuple val(tumor_normal_sample_id), path(final_svaba_somatic_sv_vcf) into svaba_sv_vcf_forDuphold
    tuple val(tumor_normal_sample_id), path(final_svaba_somatic_sv_read_support) into svaba_sv_read_support_forAnnotation

    when:
    params.svaba == "on"

    script:
    final_svaba_somatic_sv_vcf = "${tumor_normal_sample_id}.svaba.somatic.sv.vcf"
    final_svaba_somatic_sv_read_support = "${tumor_normal_sample_id}.svaba.somatic.sv.readsupp.txt"
    """
    bcftools filter \
    --output-type u \
    --exclude 'QUAL<6' \
    "${svaba_somatic_sv_vcf}" \
    | \
    bcftools filter \
    --output-type u \
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

// BCFtools Norm ~ left-align and normalize indels
process leftNormalizeSvabaVcf_bcftools {
	publishDir "${params.output_dir}/somatic/svaba", mode: 'copy', pattern: '*.{vcf.gz,tbi,txt}'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(filtered_somatic_indel_vcf), path(filtered_somatic_indel_vcf_index) from filtered_indel_vcf_forSvabaBcftools
	path ref_genome_fasta from ref_genome_fasta_file
	path ref_genome_fasta_index from ref_genome_fasta_index_file
	path ref_genome_fasta_dict from ref_genome_fasta_dict_file

	output:
	tuple val(tumor_normal_sample_id), path(final_svaba_indel_vcf) into final_svaba_indel_vcf_forConsensus
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
	--fasta-ref "${ref_genome_fasta}" \
	--output-type z \
	--output "${final_svaba_indel_vcf}" \
	"${filtered_somatic_indel_vcf}" 2>"${svaba_realign_normalize_stats}"

	tabix "${final_svaba_indel_vcf}"
	"""
}

// duphold ~ efficiently annotate SV calls with sequence depth information to reduce false positive deletion and duplication calls
process falsePostiveSvFilteringSvaba_duphold {
    publishDir "${params.output_dir}/somatic/svaba", mode: 'copy', pattern: '*.{vcf.gz}'
    tag "${tumor_normal_sample_id}"

    input:
    tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index), path(final_svaba_somatic_sv_vcf) from svaba_tumor_bam_forDuphold.join(svaba_sv_vcf_forDuphold)
    path ref_genome_fasta from ref_genome_fasta_file
    path ref_genome_fasta_index from ref_genome_fasta_index_file
    path ref_genome_fasta_dict from ref_genome_fasta_dict_file

    output:
    tuple val(tumor_normal_sample_id), path(svaba_filtered_final_sv_vcf) into svaba_filtered_final_sv_vcf_forConsensus

    when:
    params.svaba == "on"

    script:
    svaba_filtered_final_sv_vcf = "${tumor_normal_sample_id}.svaba.somatic.sv.final.vcf.gz"
    """
    duphold \
    --vcf "${final_svaba_somatic_sv_vcf}" \
    --bam "${tumor_bam}" \
    --fasta "${ref_genome_fasta}" \
    --threads ${task.cpus} \
    --output "${tumor_normal_sample_id}.svaba.somatic.sv.fpmarked.vcf"

    # Filter using recommended thresholds for DEL/DUP
    bcftools filter \
    --output-type u \
    --exclude 'INFO/SVTYPE="DEL" && FORMAT/DHFFC>0.7' \
    "${tumor_normal_sample_id}.svaba.somatic.sv.fpmarked.vcf" \
    | \
    bcftools filter \
    --output-type z \
    --exclude 'INFO/SVTYPE="DUP" && FORMAT/DHBFC<1.3' \
    --output "${svaba_filtered_final_sv_vcf}"
    """
}

// END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \\


// ~~~~~~~~~~~~~~~~ DELLY2 ~~~~~~~~~~~~~~~~~ \\
// START

// DELLY2 ~ discover structural variants using paired-ends, split-reads and read-depth
process svAndIndelCalling_delly {
	publishDir "${params.output_dir}/somatic/delly", mode: 'copy', pattern: '*.{vcf.gz,tbi}'
	tag "${tumor_normal_sample_id}"

	input:
	tuple path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index) from tumor_normal_pair_forDelly
	path ref_genome_fasta from ref_genome_fasta_file
	path ref_genome_fasta_index from ref_genome_fasta_index_file
	path ref_genome_fasta_dict from ref_genome_fasta_dict_file
	path wgs_blacklist_0based_bed from wgs_bed_blacklist_0based

	output:
	tuple val(tumor_normal_sample_id), val(tumor_id), path(delly_somatic_sv_vcf), path(delly_somatic_sv_vcf_index) into delly_sv_vcf_forPostprocessing
	tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index) into delly_tumor_bam_forDuphold

	when:
	params.delly == "on"

	script:
	tumor_id = "${tumor_bam.baseName}".replaceFirst(/\..*$/, "")
	normal_id = "${normal_bam.baseName}".replaceFirst(/\..*$/, "")
	tumor_normal_sample_id = "${tumor_id}_vs_${normal_id}"
	delly_call_parameters = params.delly_strict == "on" ? "--map-qual 30 --qual-tra 40 --mad-cutoff 15 --min-clique-size 5" : "--map-qual 1 --qual-tra 20 --mad-cutoff 9 --min-clique-size 2"
	delly_somatic_sv_vcf = "${tumor_normal_sample_id}.delly.somatic.sv.unprocessed.vcf.gz"
	delly_somatic_sv_vcf_index = "${delly_somatic_sv_vcf}.tbi"
	"""
	delly call \
	${delly_call_parameters} \
	--genome "${ref_genome_fasta}" \
	--exclude "${wgs_blacklist_0based_bed}" \
	--outfile "${tumor_normal_sample_id}.delly.sv.unfiltered.bcf" \
	"${tumor_bam}" \
	"${normal_bam}"

	touch samples.tsv
	echo "${tumor_id}\ttumor" >> samples.tsv
	echo "${normal_id}\tcontrol" >> samples.tsv

	delly filter \
	--filter somatic \
	--pass \
	--altaf 0.05 \
	--minsize 51 \
	--coverage 10 \
	--samples samples.tsv \
	"${tumor_normal_sample_id}.delly.sv.unfiltered.bcf" \
	| \
	bcftools view \
    --output-type z \
    --threads ${task.cpus} \
    --output-file "${delly_somatic_sv_vcf}"

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
    tuple val(tumor_normal_sample_id), path(final_delly_somatic_sv_vcf) into delly_sv_vcf_forDuphold
    tuple val(tumor_normal_sample_id), path(final_delly_somatic_sv_read_support) into delly_sv_read_support_forAnnotation

    when:
    params.delly == "on"

    script:
    final_delly_somatic_sv_vcf = "${tumor_normal_sample_id}.delly.somatic.sv.vcf"
    final_delly_somatic_sv_read_support = "${tumor_normal_sample_id}.delly.somatic.sv.readsupp.txt"
    """
    bcftools filter \
    --include 'INFO/MAPQ=60 || INFO/SRMAPQ=60' \
    "${delly_somatic_sv_vcf}" \
    | \
    bcftools filter \
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

// duphold ~ efficiently annotate SV calls with sequence depth information to reduce false positive deletion and duplication calls
process falsePostiveSvFilteringDelly_duphold {
    publishDir "${params.output_dir}/somatic/delly", mode: 'copy', pattern: '*.{vcf.gz}'
    tag "${tumor_normal_sample_id}"

    input:
    tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index), path(final_delly_somatic_sv_vcf) from delly_tumor_bam_forDuphold.join(delly_sv_vcf_forDuphold)
    path ref_genome_fasta from ref_genome_fasta_file
    path ref_genome_fasta_index from ref_genome_fasta_index_file
    path ref_genome_fasta_dict from ref_genome_fasta_dict_file

    output:
    tuple val(tumor_normal_sample_id), path(delly_filtered_final_sv_vcf) into delly_filtered_final_sv_vcf_forConsensus

    when:
    params.delly == "on"

    script:
    delly_filtered_final_sv_vcf = "${tumor_normal_sample_id}.delly.somatic.sv.final.vcf.gz"
    """
    duphold \
    --vcf "${final_delly_somatic_sv_vcf}" \
    --bam "${tumor_bam}" \
    --fasta "${ref_genome_fasta}" \
    --threads ${task.cpus} \
    --output "${tumor_normal_sample_id}.delly.somatic.sv.fpmarked.vcf"

    # Filter using recommended thresholds for DEL/DUP
    bcftools filter \
    --output-type u \
    --exclude 'INFO/SVTYPE="DEL" && FORMAT/DHFFC>0.7' \
    "${tumor_normal_sample_id}.delly.somatic.sv.fpmarked.vcf" \
    | \
    bcftools filter \
    --output-type z \
    --exclude 'INFO/SVTYPE="DUP" && FORMAT/DHBFC<1.3' \
    --output "${delly_filtered_final_sv_vcf}"
    """
}

// END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \\


// ~~~~~~~~~~~~~~~~ IgCaller ~~~~~~~~~~~~~~~ \\
// START

// IgCaller ~ characterize the immunoglobulin gene rearrangements and oncogenic translocations
process igRearrangementsAndTranslocations_igcaller {
    publishDir "${params.output_dir}/somatic/igcaller", mode: 'copy', pattern: '*.{tsv}'
    tag "${tumor_normal_sample_id}"

    input:
    tuple path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index) from tumor_normal_pair_forIgCaller
    path ref_genome_fasta from ref_genome_fasta_file
	path ref_genome_fasta_index from ref_genome_fasta_index_file
	path ref_genome_fasta_dict from ref_genome_fasta_dict_file
     
    output:
    tuple val(tumor_normal_sample_id), path(igcaller_oncogenic_rearrangements_tsv) into igcaller_onco_tsv_forConsensus
    path igcaller_csr_tsv
    path igcaller_igk_tsv
    path igcaller_igl_tsv
    path igcaller_igh_tsv
    path igcaller_filtered_calls_tsv

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
    --refGenome "${ref_genome_fasta}" \
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


// ~~~~~~~~~~~~~~~~ VarScan2 ~~~~~~~~~~~~~~~~ \\
// START

// VarScan somatic / SAMtools mpileup ~ heuristic/statistic approach to call SNV and indel variants
process snvAndIndelCalling_varscan {
	tag "${tumor_normal_sample_id} C=${chromosome}"

	input:
	tuple path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index) from tumor_normal_pair_forVarscanSamtoolsMpileup
	path ref_genome_fasta from ref_genome_fasta_file
	path ref_genome_fasta_index from ref_genome_fasta_index_file
	path ref_genome_fasta_dict from ref_genome_fasta_dict_file
	path wgs_bed from gatk_wgs_bed
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
	--positions "${wgs_bed}" \
	--region "${chromosome}" \
	--fasta-ref "${ref_genome_fasta}" \
	"${normal_bam}" "${tumor_bam}" \
	| \
	java -Xmx2G -XX:ParallelGCThreads=2 -jar \${VARSCAN} somatic \
	--mpileup 1 \
	--min-coverage-normal 6 \
	--min-coverage-tumor 4 \
	--min-var-freq 0.01 \
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
	java -Xmx2G -XX:ParallelGCThreads=2 -jar \${VARSCAN} processSomatic \
	"${tumor_normal_sample_id}.snv" \
	--min-tumor-freq 0.01 \
	--max-normal-freq 0.05 \
	--p-value 0.05

	bgzip < "${tumor_normal_sample_id}.snv.Somatic.hc" > "${high_confidence_snv_vcf}"
	tabix "${high_confidence_snv_vcf}"

	zcat "${raw_indel_vcf}" \
	| \
	java -Xmx2G -XX:ParallelGCThreads=2 -jar \${VARSCAN} processSomatic \
	"${tumor_normal_sample_id}.indel" \
	--min-tumor-freq 0.01 \
	--max-normal-freq 0.05 \
	--p-value 0.05

	bgzip < "${tumor_normal_sample_id}.indel.Somatic.hc" > "${high_confidence_indel_vcf}"
	tabix "${high_confidence_indel_vcf}"
	"""
}

// bam-readcount / BCFtools concat ~ generate metrics at single nucleotide positions for filtering out false positive calls
process bamReadcountForVarscanFpFilter_bamreadcount {
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(high_confidence_snv_vcf), path(high_confidence_snv_vcf_index), path(high_confidence_indel_vcf), path(high_confidence_indel_vcf_index) from bams_forVarscanBamReadcount.join(high_confidence_vcfs_forVarscanBamReadcount)
	path ref_genome_fasta from ref_genome_fasta_file
	path ref_genome_fasta_index from ref_genome_fasta_index_file
	path ref_genome_fasta_dict from ref_genome_fasta_dict_file

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
	"${ref_genome_fasta}" \
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

	java -Xmx2G -XX:ParallelGCThreads=2 -jar \$VARSCAN fpfilter \
	"${unzipped_hc_snv_vcf}" \
	"${snv_readcount_file}" \
	--min-var-count 2 \
	--min-var-freq 0.001 \
	--min-ref-basequal 25 \
	--min-var-basequal 25 \
	--output-file "${fp_filtered_snv_vcf}"

	gunzip -f "${high_confidence_indel_vcf}"

	java -Xmx2G -XX:ParallelGCThreads=2 -jar \$VARSCAN fpfilter \
	"${unzipped_hc_indel_vcf}" \
	"${indel_readcount_file}" \
	--min-var-count 2 \
	--min-var-freq 0.001 \
	--min-ref-basequal 25 \
	--min-var-basequal 25 \
	--output-file "${fp_filtered_indel_vcf}"
	"""
}

// BCFtools norm ~ split multiallelic sites into multiple rows then left-align and normalize indels
process splitMultiallelicAndLeftNormalizeVarscanVcf_bcftools {
	publishDir "${params.output_dir}/somatic/varscan", mode: 'copy', pattern: '*.{vcf.gz,tbi,txt}'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(fp_filtered_snv_vcf), path(fp_filtered_indel_vcf) from filtered_vcfs_forVarscanBcftools
	path ref_genome_fasta from ref_genome_fasta_file
	path ref_genome_fasta_index from ref_genome_fasta_index_file
	path ref_genome_fasta_dict from ref_genome_fasta_dict_file

	output:
	tuple val(tumor_normal_sample_id), path(final_varscan_snv_vcf) into final_varscan_snv_vcf_forConsensus
	tuple val(tumor_normal_sample_id), path(final_varscan_indel_vcf) into final_varscan_indel_vcf_forConsensus
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
	--multiallelics -indels \
	--output-type u \
	- 2>"${varscan_indel_multiallelics_stats}" \
	| \
	bcftools norm \
	--threads ${task.cpus} \
	--fasta-ref "${ref_genome_fasta}" \
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
	path gnomad_vcf_chromosomes1_9 from gnomad_ref_vcf_chromosomes1_9
	path gnomad_vcf_chromosomes1_9_index from gnomad_ref_vcf_chromosomes1_9_index
	path gnomad_vcf_chromosomes10_22 from gnomad_ref_vcf_chromosomes10_22
	path gnomad_vcf_chromosomes10_22_index from gnomad_ref_vcf_chromosomes10_22_index
	path gnomad_vcf_chromosomesXYM_alts from gnomad_ref_vcf_chromosomesXYM_alts
	path gnomad_vcf_chromosomesXYM_alts_index from gnomad_ref_vcf_chromosomesXYM_alts_index

	output:
	tuple path(mutect_gnomad_vcf), path(mutect_gnomad_vcf_index) into mutect_gnomad_ref_vcf_fromProcess

	when:
	params.mutect == "on" && params.mutect_ref_vcf_concatenated == "no"

	script:
	mutect_gnomad_vcf = "af-only-gnomad.hg38.vcf.gz"
	mutect_gnomad_vcf_index = "${mutect_gnomad_vcf}.tbi"
	"""
	bcftools concat \
	--threads ${task.cpus} \
	--output-type z \
	--output "${mutect_gnomad_vcf}" \
	"${gnomad_vcf_chromosomes1_9}" \
	"${gnomad_vcf_chromosomes10_22}" \
	"${gnomad_vcf_chromosomesXYM_alts}"

	tabix "${mutect_gnomad_vcf}"
	"""
}

// Depending on whether the gnomAD allele frequency reference VCF was pre-built, set the input
// channel for the for Mutect2 process
if( params.mutect_ref_vcf_concatenated == "yes" && params.mutect == "on") {
	mutect_gnomad_ref_vcf = mutect_gnomad_ref_vcf_preBuilt.combine(mutect_gnomad_ref_vcf_index_preBuilt)
}
else {
	mutect_gnomad_ref_vcf = mutect_gnomad_ref_vcf_fromProcess
}


// GATK Mutect2 ~ call somatic SNVs and indels via local assembly of haplotypes
process snvAndIndelCalling_gatk {
	tag "${tumor_normal_sample_id} C=${chromosome}"

	input:
	tuple path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(mutect_gnomad_vcf), path(mutect_gnomad_vcf_index) from tumor_normal_pair_forMutectCalling.combine(mutect_gnomad_ref_vcf)
	path ref_genome_fasta from ref_genome_fasta_file
	path ref_genome_fasta_index from ref_genome_fasta_index_file
	path ref_genome_fasta_dict from ref_genome_fasta_dict_file
	path wgs_bed from gatk_wgs_bed
	path pon_1000G from panel_of_normals_1000G
	path pon_1000G_index from panel_of_normals_1000G_index
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
	per_chromosome_bed_file = "${wgs_bed}".replaceFirst(/\.bed/, ".${chromosome}.bed")
	raw_per_chromosome_vcf = "${tumor_normal_sample_id}.${chromosome}.vcf.gz"
	raw_per_chromosome_vcf_index = "${raw_per_chromosome_vcf}.tbi"
	raw_per_chromosome_mutect_stats_file = "${tumor_normal_sample_id}.${chromosome}.vcf.gz.stats"
	"""
	grep -w '${chromosome}' "${wgs_bed}" > "${per_chromosome_bed_file}"

	gatk Mutect2 \
	--java-options "-Xmx${task.memory.toGiga()}G -Djava.io.tmpdir=." \
	--verbosity ERROR \
	--tmp-dir . \
	--native-pair-hmm-threads ${task.cpus} \
	--af-of-alleles-not-in-resource 0.00003125 \
	--seconds-between-progress-updates 600 \
	--reference "${ref_genome_fasta}" \
	--intervals "${per_chromosome_bed_file}" \
	--germline-resource "${mutect_gnomad_vcf}" \
	--panel-of-normals "${pon_1000G}" \
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

// GATK GetPileupSummaries ~ tabulates pileup metrics for inferring contamination
process pileupSummariesForMutect2Contamination_gatk {
	tag "${tumor_normal_sample_id} C=${chromosome}"

	input:
	tuple path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index) from tumor_normal_pair_forMutectPileup
	path wgs_bed from gatk_wgs_bed
	path exac_common_sites_vcf from exac_common_sites_ref_vcf
	path exac_common_sites_vcf_index from exac_common_sites_ref_vcf_index
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
	per_chromosome_bed_file = "${wgs_bed}".replaceFirst(/\.bed/, ".${chromosome}.bed")
	per_chromosome_tumor_pileup = "${tumor_id}.${chromosome}.pileup"
	per_chromosome_normal_pileup = "${normal_id}.${chromosome}.pileup"
	"""
	grep -w '${chromosome}' "${wgs_bed}" > "${per_chromosome_bed_file}"

	gatk GetPileupSummaries \
	--java-options "-Xmx${task.memory.toGiga()}G -Djava.io.tmpdir=." \
	--verbosity ERROR \
	--tmp-dir . \
	--intervals "${per_chromosome_bed_file}" \
	--variant "${exac_common_sites_vcf}" \
	--input "${tumor_bam}" \
	--output "${per_chromosome_tumor_pileup}"

	gatk GetPileupSummaries \
	--java-options "-Xmx${task.memory.toGiga()}G -Djava.io.tmpdir=." \
	--verbosity ERROR \
	--tmp-dir . \
	--intervals "${per_chromosome_bed_file}" \
	--variant "${exac_common_sites_vcf}" \
	--input "${normal_bam}" \
	--output "${per_chromosome_normal_pileup}"
	"""
}

// GATK GatherPileupSummaries ~ combine tumor pileup tables for inferring contamination
process gatherTumorPileupSummariesForMutect2Contamination_gatk {
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(per_chromosome_tumor_pileup) from per_chromosome_tumor_pileups_forMutectPileupGather.groupTuple()
	path ref_genome_fasta_dict from ref_genome_fasta_dict_file

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
	--sequence-dictionary "${ref_genome_fasta_dict}" \
	${per_chromosome_tumor_pileup.collect { "--I $it " }.join()} \
	--O "${tumor_pileup}"
	"""
}

// GATK GatherPileupSummaries ~ combine normal pileup tables for inferring contamination
process gatherNormalPileupSummariesForMutect2Contamination_gatk {
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(per_chromosome_normal_pileup) from per_chromosome_normal_pileups_forMutectPileupGather.groupTuple()
	path ref_genome_fasta_dict from ref_genome_fasta_dict_file

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
	--sequence-dictionary "${ref_genome_fasta_dict}" \
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

// GATK FilterMutectCalls ~ filter somatic SNVs and indels called by Mutect2
process mutect2VariantFiltration_gatk {
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(merged_raw_vcf), path(merged_raw_vcf_index), path(merged_mutect_stats_file), path(mutect_contamination_file) from merged_raw_vcfs_forMutectFilter.join( merged_mutect_stats_file_forMutectFilter ).join( contamination_file_forMutectFilter )
	path ref_genome_fasta from ref_genome_fasta_file
	path ref_genome_fasta_index from ref_genome_fasta_index_file
	path ref_genome_fasta_dict from ref_genome_fasta_dict_file

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
	--unique-alt-read-count 4 \
	--reference "${ref_genome_fasta}" \
	--stats "${merged_mutect_stats_file}" \
	--variant "${merged_raw_vcf}" \
	--contamination-table "${mutect_contamination_file}" \
	--output "${filtered_vcf}" \
	--filtering-stats "${filter_stats_file}"
	"""
}

// BCFtools Norm ~ split multiallelic sites into multiple rows then left-align and normalize indels
process splitMultiallelicAndLeftNormalizeMutect2Vcf_bcftools {
	publishDir "${params.output_dir}/somatic/mutect", mode: 'copy', pattern: '*.{txt}'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(filtered_vcf), path(filtered_vcf_index) from filtered_vcf_forMutectBcftools
	path ref_genome_fasta from ref_genome_fasta_file
	path ref_genome_fasta_index from ref_genome_fasta_index_file
	path ref_genome_fasta_dict from ref_genome_fasta_dict_file

	output:
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
	--output-type u \
	- 2>"${mutect_multiallelics_stats}" \
	| \
	bcftools norm \
	--threads ${task.cpus} \
	--fasta-ref "${ref_genome_fasta}" \
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
	tuple val(tumor_normal_sample_id), path(final_mutect_snv_vcf) into final_mutect_snv_vcf_forConsensus
	tuple val(tumor_normal_sample_id), path(final_mutect_indel_vcf) into final_mutect_indel_vcf_forConsensus

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


// ~~~~~~~~~~~~~~~~ Conpair ~~~~~~~~~~~~~~~~ \\
// START

// Conpair run_gatk_pileup_for_sample ~ generate GATK pileups the tumor and normal BAMs separately
process bamPileupForConpair_conpair {
	tag "${tumor_normal_sample_id}"

	input:
	tuple path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index) from tumor_normal_pair_forConpairPileup
	path ref_genome_fasta from ref_genome_fasta_file
	path ref_genome_fasta_index from ref_genome_fasta_index_file
	path ref_genome_fasta_dict from ref_genome_fasta_dict_file

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
	--reference "${ref_genome_fasta}" \
	--markers \${CONPAIR_DIR}"${hg38_ref_genome_markers}.bed"

	\${CONPAIR_DIR}/scripts/run_gatk_pileup_for_sample.py \
	--xmx_jav "${task.memory.toGiga()}g" \
	--bam "${normal_bam}" \
	--outfile "${normal_pileup}" \
	--reference "${ref_genome_fasta}" \
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
	tuple val(tumor_normal_sample_id), path(conpair_contamination_file) into normal_contamination_forCaveman

	when:
	params.conpair == "on"
	
	script:
	conpair_concordance_file = "${tumor_normal_sample_id}.conpair.concordance.txt"
	conpair_contamination_file = "${tumor_normal_sample_id}.conpair.contamination.txt"
	hg38_ref_genome_markers = "/data/markers/GRCh38.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.liftover"
	"""
	\${CONPAIR_DIR}/scripts/verify_concordance.py \
	--min_cov ${params.conpair_min_cov} \
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
	params.caveman == "on" && params.battenberg == "on" && params.manta == "on" && params.conpair == "on"

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

// CaVEMan setup ~ generate configuration files for all execution steps
process setup_caveman {
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(tumor_cnv_profile), path(germline_indel_bed), path(germline_indel_bed_index), path(normal_contamination_file) from bams_forCaveman.join(tumor_cnv_profile_forCaveman).join(germline_indel_bed_forCaveman).join(normal_contamination_forCaveman)
	path ref_genome_fasta from ref_genome_fasta_file
	path ref_genome_fasta_index from ref_genome_fasta_index_file
	path ref_genome_fasta_dict from ref_genome_fasta_dict_file
	path blacklist from wes_bed_blacklist_0based
	path unmatched_normal from unmatched_normal_bed
	path unmatched_normal_index from unmatched_normal_bed_index
	path centromeric_repeats from centromeric_repeats_bed
	path centromeric_repeats_index from centromeric_repeats_bed_index
	path simple_repeats from simple_repeats_bed
	path simple_repeats_index from simple_repeats_bed_index
	path dbsnp from dbsnp_bed
	path dbsnp_index from dbsnp_bed_index

	output:
	tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(tumor_cnv_profile_bed), path(normal_cnv_profile_bed), path(germline_indel_bed), path(germline_indel_bed_index), path(normal_contamination_file), path(postprocessing_config_file), path(config_file), path(alg_bean_file) into setup_forCavemanSplit, setup_forCavemanConcat, setup_forCavemanMstep, setup_forCavemanMerge, setup_forCavemanEstep, setup_forCavemanFlag, setup_forCavemanResults

	when:
	params.caveman == "on" && params.battenberg == "on" && params.manta == "on" && params.conpair == "on"

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
	echo "vcfUnmatchedMinSamplePct=1.000" >> "${postprocessing_config_file}"
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
	echo "centromericRepeatBed=${centromeric_repeats}" >> "${postprocessing_config_file}"
	echo "simpleRepeatBed=${simple_repeats}" >> "${postprocessing_config_file}"
	echo "snpBed=${dbsnp}" >> "${postprocessing_config_file}"
	echo "" >> "${postprocessing_config_file}"

	grep -v 'startpos' "${tumor_cnv_profile}" \
	| \
	awk 'BEGIN {OFS="\t"} {print \$1,\$2,\$3,\$8+\$9}' > "${tumor_cnv_profile_bed}"

	touch "${normal_cnv_profile_bed}"

	normal_contamination=\$(awk -v contam_percent=\$(grep 'Tumor' ${normal_contamination_file} | cut -d ' ' -f 5 | sed 's|%||') -v denom=100 'BEGIN {print  (contam_percent / denom)}')

	caveman.pl \
	-outdir . \
	-reference "${ref_genome_fasta_index}" \
	-ignore-file "${blacklist}" \
	-tumour-bam "${tumor_bam}" \
	-normal-bam "${normal_bam}" \
	-tumour-cn "${tumor_cnv_profile_bed}" \
	-normal-cn "${normal_cnv_profile_bed}" \
	-norm-cn-default 2 \
	-species Homo_sapiens \
	-species-assembly GRCh38 \
	-flag-bed-files . \
	-germline-indel "${germline_indel_bed}" \
	-unmatched-vcf "${unmatched_normal}" \
	-seqType genome \
	-threads ${task.cpus} \
	-normal-contamination \${normal_contamination} \
	-flagConfig "${postprocessing_config_file}" \
	-read-count 3500000 \
	-process setup \
	-index 1
	"""
}

// CaVEMan split ~ split the genome into chunks by readsize and hard stop forced by contig ends
process split_caveman {
	tag "${tumor_normal_sample_id} C=${chromosome}"

	input:
	tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(tumor_cnv_profile_bed), path(normal_cnv_profile_bed), path(germline_indel_bed), path(germline_indel_bed_index), path(normal_contamination_file), path(postprocessing_config_file), path(config_file), path(alg_bean_file) from setup_forCavemanSplit
	path ref_genome_fasta from ref_genome_fasta_file
	path ref_genome_fasta_index from ref_genome_fasta_index_file
	path ref_genome_fasta_dict from ref_genome_fasta_dict_file
	path blacklist from wes_bed_blacklist_0based
	path unmatched_normal from unmatched_normal_bed
	path unmatched_normal_index from unmatched_normal_bed_index
	path centromeric_repeats from centromeric_repeats_bed
	path centromeric_repeats_index from centromeric_repeats_bed_index
	path simple_repeats from simple_repeats_bed
	path simple_repeats_index from simple_repeats_bed_index
	path dbsnp from dbsnp_bed
	path dbsnp_index from dbsnp_bed_index
	each chromosome from chromosome_list_forCavemanSplit

	output:
	tuple val(tumor_normal_sample_id), path(split_list_per_chromosome), path(read_position_per_chromosome) into split_per_chromosome_forCavemanConcat, split_per_chromosome_forCavemanMstep, split_per_chromosome_forCavemanMerge, split_per_chromosome_forCavemanEstep, split_per_chromosome_forCavemanFlag, split_per_chromosome_forCavemanResults

	when:
	params.caveman == "on" && params.battenberg == "on" && params.manta == "on" && params.conpair == "on"

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

	normal_contamination=\$(awk -v contam_percent=\$(grep 'Tumor' ${normal_contamination_file} | cut -d ' ' -f 5 | sed 's|%||') -v denom=100 'BEGIN {print  (contam_percent / denom)}')

	i=\$(grep -wn "${chromosome}" "${ref_genome_fasta_index}" | cut -f 1 | cut -d ':' -f 1)

	caveman.pl \
	-outdir . \
	-reference "${ref_genome_fasta_index}" \
	-ignore-file "${blacklist}" \
	-tumour-bam "${tumor_bam}" \
	-normal-bam "${normal_bam}" \
	-tumour-cn "${tumor_cnv_profile_bed}" \
	-normal-cn "${normal_cnv_profile_bed}" \
	-norm-cn-default 2 \
	-species Homo_sapiens \
	-species-assembly GRCh38 \
	-flag-bed-files . \
	-germline-indel "${germline_indel_bed}" \
	-unmatched-vcf "${unmatched_normal}" \
	-seqType genome \
	-threads ${task.cpus} \
	-normal-contamination \${normal_contamination} \
	-flagConfig "${postprocessing_config_file}" \
	-read-count 3500000 \
	-process split \
	-index \${i}
	"""
}

// CaVEMan split_concat ~ concatenate the split file sections into a single split section reference file
process splitConcat_caveman {
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(tumor_cnv_profile_bed), path(normal_cnv_profile_bed), path(germline_indel_bed), path(germline_indel_bed_index), path(normal_contamination_file), path(postprocessing_config_file), path(config_file), path(alg_bean_file), path(split_list_per_chromosome), path(read_position_per_chromosome) from setup_forCavemanConcat.join(split_per_chromosome_forCavemanConcat.groupTuple())
	path ref_genome_fasta from ref_genome_fasta_file
	path ref_genome_fasta_index from ref_genome_fasta_index_file
	path ref_genome_fasta_dict from ref_genome_fasta_dict_file
	path blacklist from wes_bed_blacklist_0based
	path unmatched_normal from unmatched_normal_bed
	path unmatched_normal_index from unmatched_normal_bed_index
	path centromeric_repeats from centromeric_repeats_bed
	path centromeric_repeats_index from centromeric_repeats_bed_index
	path simple_repeats from simple_repeats_bed
	path simple_repeats_index from simple_repeats_bed_index
	path dbsnp from dbsnp_bed
	path dbsnp_index from dbsnp_bed_index

	output:
	tuple val(tumor_normal_sample_id), path(split_list), path(index_step_list) into split_concat_forCavemanMstep, split_concat_forCavemanMerge, split_concat_forCavemanEstep, split_concat_forCavemanFlag, split_concat_forCavemanResults

	when:
	params.caveman == "on" && params.battenberg == "on" && params.manta == "on" && params.conpair == "on"

	script:
	split_list = "tmpCaveman/splitList"
	index_step_list = "${tumor_normal_sample_id}.step_list"
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
	mv splitList.chr* tmpCaveman/
	mv readpos.chr* tmpCaveman/

	normal_contamination=\$(awk -v contam_percent=\$(grep 'Tumor' ${normal_contamination_file} | cut -d ' ' -f 5 | sed 's|%||') -v denom=100 'BEGIN {print  (contam_percent / denom)}')

	caveman.pl \
	-outdir . \
	-reference "${ref_genome_fasta_index}" \
	-ignore-file "${blacklist}" \
	-tumour-bam "${tumor_bam}" \
	-normal-bam "${normal_bam}" \
	-tumour-cn "${tumor_cnv_profile_bed}" \
	-normal-cn "${normal_cnv_profile_bed}" \
	-norm-cn-default 2 \
	-species Homo_sapiens \
	-species-assembly GRCh38 \
	-flag-bed-files . \
	-germline-indel "${germline_indel_bed}" \
	-unmatched-vcf "${unmatched_normal}" \
	-seqType genome \
	-threads ${task.cpus} \
	-normal-contamination \${normal_contamination} \
	-flagConfig "${postprocessing_config_file}" \
	-read-count 3500000 \
	-process split_concat \
	-index 1

	seq 1 \$(cat ${split_list} | wc -l) > "${index_step_list}"
	"""
}

setup_forCavemanMstep.join(split_per_chromosome_forCavemanMstep.groupTuple())
	.join(split_concat_forCavemanMstep)
	.splitText(elem: 19, by:1)
	.set{ input_forCavemanMstep }

// CaVEMan mstep ~ build a profile of each split section of the genome using various covariates
process mstep_caveman {
	tag "${tumor_normal_sample_id} IDX=${index}"

	input:
	tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(tumor_cnv_profile_bed), path(normal_cnv_profile_bed), path(germline_indel_bed), path(germline_indel_bed_index), path(normal_contamination_file), path(postprocessing_config_file), path(config_file), path(alg_bean_file), path(split_list_per_chromosome), path(read_position_per_chromosome), path(split_list), val(index_dirty) from input_forCavemanMstep
	path ref_genome_fasta from ref_genome_fasta_file
	path ref_genome_fasta_index from ref_genome_fasta_index_file
	path ref_genome_fasta_dict from ref_genome_fasta_dict_file
	path blacklist from wes_bed_blacklist_0based
	path unmatched_normal from unmatched_normal_bed
	path unmatched_normal_index from unmatched_normal_bed_index
	path centromeric_repeats from centromeric_repeats_bed
	path centromeric_repeats_index from centromeric_repeats_bed_index
	path simple_repeats from simple_repeats_bed
	path simple_repeats_index from simple_repeats_bed_index
	path dbsnp from dbsnp_bed
	path dbsnp_index from dbsnp_bed_index

	output:
	tuple val(tumor_normal_sample_id), path(mstep_results_directory_per_index) into mstep_covs_forCavemanMerge, mstep_covs_forCavemanEstep

	when:
	params.caveman == "on" && params.battenberg == "on" && params.manta == "on" && params.conpair == "on"

	script:
	index = "${index_dirty}".replaceFirst("\\s","")
	mstep_results_directory_per_index = "results_mstep_${index}"
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
	mv splitList.chr* tmpCaveman/
	mv readpos.chr* tmpCaveman/
	mv "${split_list}" tmpCaveman/

	normal_contamination=\$(awk -v contam_percent=\$(grep 'Tumor' ${normal_contamination_file} | cut -d ' ' -f 5 | sed 's|%||') -v denom=100 'BEGIN {print  (contam_percent / denom)}')

	caveman.pl \
	-outdir . \
	-reference "${ref_genome_fasta_index}" \
	-ignore-file "${blacklist}" \
	-tumour-bam "${tumor_bam}" \
	-normal-bam "${normal_bam}" \
	-tumour-cn "${tumor_cnv_profile_bed}" \
	-normal-cn "${normal_cnv_profile_bed}" \
	-norm-cn-default 2 \
	-species Homo_sapiens \
	-species-assembly GRCh38 \
	-flag-bed-files . \
	-germline-indel "${germline_indel_bed}" \
	-unmatched-vcf "${unmatched_normal}" \
	-seqType genome \
	-threads ${task.cpus} \
	-normal-contamination \${normal_contamination} \
	-flagConfig "${postprocessing_config_file}" \
	-read-count 3500000 \
	-process mstep \
	-index "${index}"

	mkdir "${mstep_results_directory_per_index}"
	cp -a tmpCaveman/results/* "${mstep_results_directory_per_index}"
	"""
}

// CaVEMan merge ~ build a single profile of the whole genome using the mstep covariates
process merge_caveman {
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(tumor_cnv_profile_bed), path(normal_cnv_profile_bed), path(germline_indel_bed), path(germline_indel_bed_index), path(normal_contamination_file), path(postprocessing_config_file), path(config_file), path(alg_bean_file), path(split_list_per_chromosome), path(read_position_per_chromosome), path(split_list), path(index_step_list), path(mstep_results_directory_per_index) from setup_forCavemanMerge.join(split_per_chromosome_forCavemanMerge.groupTuple()).join(split_concat_forCavemanMerge).join(mstep_covs_forCavemanMerge.groupTuple())
	path ref_genome_fasta from ref_genome_fasta_file
	path ref_genome_fasta_index from ref_genome_fasta_index_file
	path ref_genome_fasta_dict from ref_genome_fasta_dict_file
	path blacklist from wes_bed_blacklist_0based
	path unmatched_normal from unmatched_normal_bed
	path unmatched_normal_index from unmatched_normal_bed_index
	path centromeric_repeats from centromeric_repeats_bed
	path centromeric_repeats_index from centromeric_repeats_bed_index
	path simple_repeats from simple_repeats_bed
	path simple_repeats_index from simple_repeats_bed_index
	path dbsnp from dbsnp_bed
	path dbsnp_index from dbsnp_bed_index

	output:
	tuple val(tumor_normal_sample_id), path(covariate_file), path(probabilities_file) into merged_covs_probs_forCavemanEstep

	when:
	params.caveman == "on" && params.battenberg == "on" && params.manta == "on" && params.conpair == "on"

	script:
	covariate_file = "tmpCaveman/cov_arr"
	probabilities_file = "tmpCaveman/prob_arr"
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
	mv splitList.chr* tmpCaveman/
	mv readpos.chr* tmpCaveman/
	mv "${split_list}" tmpCaveman/
	mkdir -p tmpCaveman/results
	cp -a results_mstep_*/* tmpCaveman/results/

	normal_contamination=\$(awk -v contam_percent=\$(grep 'Tumor' ${normal_contamination_file} | cut -d ' ' -f 5 | sed 's|%||') -v denom=100 'BEGIN {print  (contam_percent / denom)}')

	caveman.pl \
	-outdir . \
	-reference "${ref_genome_fasta_index}" \
	-ignore-file "${blacklist}" \
	-tumour-bam "${tumor_bam}" \
	-normal-bam "${normal_bam}" \
	-tumour-cn "${tumor_cnv_profile_bed}" \
	-normal-cn "${normal_cnv_profile_bed}" \
	-norm-cn-default 2 \
	-species Homo_sapiens \
	-species-assembly GRCh38 \
	-flag-bed-files . \
	-germline-indel "${germline_indel_bed}" \
	-unmatched-vcf "${unmatched_normal}" \
	-seqType genome \
	-threads ${task.cpus} \
	-normal-contamination \${normal_contamination} \
	-flagConfig "${postprocessing_config_file}" \
	-read-count 3500000 \
	-process merge \
	-index 1
	"""
}

setup_forCavemanEstep.join(split_per_chromosome_forCavemanEstep.groupTuple())
	.join(mstep_covs_forCavemanEstep.groupTuple())
	.join(merged_covs_probs_forCavemanEstep)
	.join(split_concat_forCavemanEstep)
	.splitText(elem: 22, by:1)
	.set{ input_forCavemanEstep }

// CaVEMan estep ~ assign probability to genotype at each position using mstep profile, sequence data, and copy number information
process snvCalling_caveman {
	tag "${tumor_normal_sample_id} IDX=${index}"

	input:
	tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(tumor_cnv_profile_bed), path(normal_cnv_profile_bed), path(germline_indel_bed), path(germline_indel_bed_index), path(normal_contamination_file), path(postprocessing_config_file), path(config_file), path(alg_bean_file), path(split_list_per_chromosome), path(read_position_per_chromosome), path(mstep_results_directory_per_index), path(covariate_file), path(probabilities_file), path(split_list), val(index_dirty) from input_forCavemanEstep
	path ref_genome_fasta from ref_genome_fasta_file
	path ref_genome_fasta_index from ref_genome_fasta_index_file
	path ref_genome_fasta_dict from ref_genome_fasta_dict_file
	path blacklist from wes_bed_blacklist_0based
	path unmatched_normal from unmatched_normal_bed
	path unmatched_normal_index from unmatched_normal_bed_index
	path centromeric_repeats from centromeric_repeats_bed
	path centromeric_repeats_index from centromeric_repeats_bed_index
	path simple_repeats from simple_repeats_bed
	path simple_repeats_index from simple_repeats_bed_index
	path dbsnp from dbsnp_bed
	path dbsnp_index from dbsnp_bed_index

	output:
	tuple val(tumor_normal_sample_id), path(estep_results_directory_per_index) into raw_vcfs_forCavemanFlag

	when:
	params.caveman == "on" && params.battenberg == "on" && params.manta == "on" && params.conpair == "on"

	script:
	index = "${index_dirty}".replaceFirst("\\s","")
	estep_results_directory_per_index = "results_estep_${index}"
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
	mv splitList.chr* tmpCaveman/
	mv readpos.chr* tmpCaveman/
	mv "${split_list}" tmpCaveman/
	mv "${covariate_file}" tmpCaveman/
	mv "${probabilities_file}" tmpCaveman/
	mkdir -p tmpCaveman/results
	cp -a results_mstep_*/* tmpCaveman/results/

	normal_contamination=\$(awk -v contam_percent=\$(grep 'Tumor' ${normal_contamination_file} | cut -d ' ' -f 5 | sed 's|%||') -v denom=100 'BEGIN {print  (contam_percent / denom)}')

	caveman.pl \
	-outdir . \
	-reference "${ref_genome_fasta_index}" \
	-ignore-file "${blacklist}" \
	-tumour-bam "${tumor_bam}" \
	-normal-bam "${normal_bam}" \
	-tumour-cn "${tumor_cnv_profile_bed}" \
	-normal-cn "${normal_cnv_profile_bed}" \
	-norm-cn-default 2 \
	-species Homo_sapiens \
	-species-assembly GRCh38 \
	-flag-bed-files . \
	-germline-indel "${germline_indel_bed}" \
	-unmatched-vcf "${unmatched_normal}" \
	-seqType genome \
	-threads ${task.cpus} \
	-normal-contamination \${normal_contamination} \
	-flagConfig "${postprocessing_config_file}" \
	-read-count 3500000 \
	-process estep \
	-index "${index}"

	mkdir "${estep_results_directory_per_index}"
	mv tmpCaveman/results "${estep_results_directory_per_index}"
	"""
}

setup_forCavemanFlag.join(split_per_chromosome_forCavemanFlag.groupTuple())
	.join(raw_vcfs_forCavemanFlag.groupTuple())
	.join(split_concat_forCavemanFlag)
	.splitText(elem: 20, by:1)
	.set{ input_forCavemanFlag }

// CaVEMan flag ~ apply filtering on raw VCF calls generated using CaVEMan
process flag_caveman {
	tag "${tumor_normal_sample_id} IDX=${index}"

	input:
	tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(tumor_cnv_profile_bed), path(normal_cnv_profile_bed), path(germline_indel_bed), path(germline_indel_bed_index), path(normal_contamination_file), path(postprocessing_config_file), path(config_file), path(alg_bean_file), path(split_list_per_chromosome), path(read_position_per_chromosome), path(estep_results_directory_per_index), path(split_list), val(index_dirty) from input_forCavemanFlag
	path ref_genome_fasta from ref_genome_fasta_file
	path ref_genome_fasta_index from ref_genome_fasta_index_file
	path ref_genome_fasta_dict from ref_genome_fasta_dict_file
	path blacklist from wes_bed_blacklist_0based
	path unmatched_normal from unmatched_normal_bed
	path unmatched_normal_index from unmatched_normal_bed_index
	path centromeric_repeats from centromeric_repeats_bed
	path centromeric_repeats_index from centromeric_repeats_bed_index
	path simple_repeats from simple_repeats_bed
	path simple_repeats_index from simple_repeats_bed_index
	path dbsnp from dbsnp_bed
	path dbsnp_index from dbsnp_bed_index

	output:
	tuple val(tumor_normal_sample_id), path(flag_results_directory_per_index) into postprocessing_output_forCavemanResults

	when:
	params.caveman == "on" && params.battenberg == "on" && params.manta == "on" && params.conpair == "on"

	script:
	index = "${index_dirty}".replaceFirst("\\s","")
	flag_results_directory_per_index = "results_flag_${index}"
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
	mv splitList.chr* tmpCaveman/
	mv readpos.chr* tmpCaveman/
	mv "${split_list}" tmpCaveman/

	normal_contamination=\$(awk -v contam_percent=\$(grep 'Tumor' ${normal_contamination_file} | cut -d ' ' -f 5 | sed 's|%||') -v denom=100 'BEGIN {print  (contam_percent / denom)}')

	for input_vcf in `ls -1v results_estep_${index}/*/*/*.muts.vcf.gz`;
		do
			flagged_vcf=\$(echo \$input_vcf | sed 's|.muts.vcf.gz|.flagged.muts.vcf|')

			cgpFlagCaVEMan.pl \
			--input \${input_vcf} \
			--outFile \${flagged_vcf} \
			--species Homo_sapiens \
			--species-assembly GRCh38 \
			--tumBam "${tumor_bam}" \
			--normBam "${normal_bam}" \
			--bedFileLoc . \
			--indelBed "${germline_indel_bed}" \
			--unmatchedVCFLoc . \
			--reference "${ref_genome_fasta_index}" \
			--flagConfig "${postprocessing_config_file}" \
			--studyType WGS

			corrected_header_vcf=\$(echo \${flagged_vcf} | sed 's|.muts.vcf|.fixedheader.muts.vcf|')

			grep -v '##vcfProcessLog' \${flagged_vcf} > \${corrected_header_vcf}
		done

	mkdir -p "${flag_results_directory_per_index}"
	cp -a results_estep_${index}/* "${flag_results_directory_per_index}"
	"""
}

// CaVEMan merge_results ~ merge all per split index VCFs into a single final result
process mergeResults_caveman {
	publishDir "${params.output_dir}/somatic/caveman", mode: 'copy', pattern: '*.{vcf.gz,tbi}'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(tumor_cnv_profile_bed), path(normal_cnv_profile_bed), path(germline_indel_bed), path(germline_indel_bed_index), path(normal_contamination_file), path(postprocessing_config_file), path(config_file), path(alg_bean_file), path(split_list_per_chromosome), path(read_position_per_chromosome), path(split_list), path(index_step_list), path(flag_results_directory_per_index) from setup_forCavemanResults.join(split_per_chromosome_forCavemanResults.groupTuple()).join(split_concat_forCavemanResults).join(postprocessing_output_forCavemanResults.groupTuple())
	path ref_genome_fasta from ref_genome_fasta_file
	path ref_genome_fasta_index from ref_genome_fasta_index_file
	path ref_genome_fasta_dict from ref_genome_fasta_dict_file
	path blacklist from wes_bed_blacklist_0based
	path unmatched_normal from unmatched_normal_bed
	path unmatched_normal_index from unmatched_normal_bed_index
	path centromeric_repeats from centromeric_repeats_bed
	path centromeric_repeats_index from centromeric_repeats_bed_index
	path simple_repeats from simple_repeats_bed
	path simple_repeats_index from simple_repeats_bed_index
	path dbsnp from dbsnp_bed
	path dbsnp_index from dbsnp_bed_index

	output:
	tuple val(tumor_normal_sample_id), path(final_caveman_snv_vcf) into final_caveman_snv_vcf_forConsensus

	when:
	params.caveman == "on" && params.battenberg == "on" && params.manta == "on" && params.conpair == "on"

	script:
	final_caveman_snv_vcf = "${tumor_normal_sample_id}.caveman.somatic.snv.vcf.gz"
	final_caveman_snv_vcf_index = "${final_caveman_snv_vcf}.tbi"
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
	mv splitList.chr* tmpCaveman/
	mv readpos.chr* tmpCaveman/
	mkdir -p tmpCaveman/results
	cp -a results_flag_*/*/* tmpCaveman/results/

	mergeCavemanResults \
	--output "${tumor_normal_sample_id}.caveman.somatic.vcf" \
	--splitlist "${split_list}" \
	-f tmpCaveman/results/%/%.flagged.fixedheader.muts.vcf

	grep -E '^#' "${tumor_normal_sample_id}.caveman.somatic.vcf" \
	| \
	bgzip > "${final_caveman_snv_vcf}"

	tabix "${final_caveman_snv_vcf}"
	"""
}

// END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \\


// ~~~~~~ UNION CONSENSUS SNV/INDEL ~~~~~~~~ \\
// START

// devgru ~ merge VCF files by calls, generating a union followed by a consensus pass
process unionAndConsensusSnvCalls_devgru {
    publishDir "${params.output_dir}/somatic/consensus/${tumor_normal_sample_id}", mode: 'copy'
    tag "${tumor_normal_sample_id}"

    input:
    tuple val(tumor_normal_sample_id), path(final_varscan_snv_vcf), path(final_mutect_snv_vcf), path(final_strelka_snv_vcf) from final_varscan_snv_vcf_forConsensus.join(final_mutect_snv_vcf_forConsensus).join(final_strelka_snv_vcf_forConsensus)
    path gene_gtf from ensembl_gtf_file

    output:
    path hq_union_consensus_snv_table

    when:
    params.varscan == "on" && params.mutect == "on" && params.strelka == "on"

    script:
    hq_union_consensus_snv_table = "${tumor_normal_sample_id}.hq.union.consensus.somatic.snv.txt.gz"
    """
    mkdir -p results/

    Rscript --vanilla ${workflow.projectDir}/bin/snvindel_union_consensus_polisher.R \
    ./ \
    snv \
    results/ \
    "${gene_gtf}" \
    ${task.cpus}

    mv results/* .
    """
}

// devgru ~ merge VCF files by calls, generating a union followed by a consensus pass
process unionAndConsensusIndelCalls_devgru {
    publishDir "${params.output_dir}/somatic/consensus/${tumor_normal_sample_id}", mode: 'copy'
    tag "${tumor_normal_sample_id}"

    input:
    tuple val(tumor_normal_sample_id), path(final_varscan_indel_vcf), path(final_mutect_indel_vcf), path(final_strelka_indel_vcf), path(final_svaba_indel_vcf) from final_varscan_indel_vcf_forConsensus.join(final_mutect_indel_vcf_forConsensus).join(final_strelka_indel_vcf_forConsensus).join(final_svaba_indel_vcf_forConsensus)
    path gene_gtf from ensembl_gtf_file

    output:
    path hq_union_consensus_indel_table

    when:
    params.varscan == "on" && params.mutect == "on" && params.strelka == "on" && params.svaba == "on"

    script:
    hq_union_consensus_indel_table = "${tumor_normal_sample_id}.hq.union.consensus.somatic.indel.txt.gz"
    """
    mkdir -p results/

    Rscript --vanilla ${workflow.projectDir}/bin/snvindel_union_consensus_polisher.R \
    ./ \
    indel \
    results/ \
    "${gene_gtf}" \
    ${task.cpus}

    mv results/* .
    """
}

// END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \\


// ~~~~~~~ UNION CONSENSUS CNV BED ~~~~~~~~~ \\
// START

// BEDtools unionbedg 2-way ~ transform CNV output into BED files then generate merged CNV segment file
process twoWayMergeAndGenerateConsensusCnvCalls_bedtools {
    publishDir "${params.output_dir}/somatic/consensus/${tumor_normal_sample_id}", mode: 'copy', pattern: '*.{merged.bed}'
    tag "${tumor_normal_sample_id}"

    input:
    tuple val(tumor_normal_sample_id), path(battenberg_somatic_cnv_bed), path(battenberg_somatic_alleles_bed), path(facets_somatic_cnv_bed), path(facets_somatic_alleles_bed) from final_battenberg_cnv_profile_forConsensus.join(final_facets_cnv_profile_forConsensus)

    output:
    tuple val(tumor_normal_sample_id), path(two_way_consensus_merged_cnv_alleles_bed) into two_way_consensus_cnv_and_allele_bed_forConsensusCnvTransform

    when:
    params.battenberg == "on" && params.facets == "on"

    script:
    two_way_merged_cnv_bed = "${tumor_normal_sample_id}.merged.somatic.cnv.bed"
    two_way_consensus_cnv_bed = "${tumor_normal_sample_id}.consensus.somatic.cnv.bed"
    two_way_merged_alleles_bed = "${tumor_normal_sample_id}.merged.somatic.alleles.bed"
    two_way_consensus_alleles_bed = "${tumor_normal_sample_id}.consensus.somatic.alleles.bed"
    two_way_consensus_merged_cnv_alleles_bed = "${tumor_normal_sample_id}.hq.union.consensus.somatic.cnv.bed"
    """
    ### Create consensus total copy number file ###
    bedtools unionbedg \
    -filler . \
    -i "${battenberg_somatic_cnv_bed}" "${facets_somatic_cnv_bed}" \
    -header \
    -names battenberg_total_cn facets_total_cn > "${two_way_merged_cnv_bed}"

    #two_way_consensus_cnv_generator.py \
    #<(grep -v 'chrom' "${two_way_merged_cnv_bed}") \
    #"${two_way_consensus_cnv_bed}"

    ### Create consensus major and minor allele file ###
    bedtools unionbedg \
    -filler . \
    -i "${battenberg_somatic_alleles_bed}" "${facets_somatic_alleles_bed}" \
    -header \
    -names battenberg_major_minor_alleles facets_major_minor_alleles > "${two_way_merged_alleles_bed}"

    #two_way_consensus_allele_generator.py \
    #<(grep -v 'chrom' "${two_way_merged_alleles_bed}") \
    #"${two_way_consensus_alleles_bed}"

    ### Merge both consensus CNV and called alleles per segment ###
    paste "${two_way_merged_cnv_bed}" <(cut -f 4-8 "${two_way_merged_alleles_bed}") \
    | \
    awk 'BEGIN {OFS="\t"} {print \$1,\$2,\$3,\$4,\$6,\$5,\$7}' > "${two_way_consensus_merged_cnv_alleles_bed}"
    """
}

// END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \\


// ~~~~~~~~ UNION CONSENSUS SV VCF ~~~~~~~~~ \\
// START

// gGnome ~ merge SV VCF files to generate a consensus
process mergeAndGenerateConsensusSvCalls_ggnome {
    publishDir "${params.output_dir}/somatic/consensus/${tumor_normal_sample_id}", mode: 'copy', pattern: '*.{bedpe,pdf}'
    tag  "${tumor_normal_sample_id}"

    input:
    tuple val(tumor_normal_sample_id), path(manta_filtered_final_sv_vcf), path(svaba_filtered_final_sv_vcf), path(delly_filtered_final_sv_vcf), path(igcaller_oncogenic_rearrangements_tsv) from manta_filtered_final_sv_vcf_forConsensus.join(svaba_filtered_final_sv_vcf_forConsensus).join(delly_filtered_final_sv_vcf_forConsensus).join(igcaller_onco_tsv_forConsensus)

    output:
    path hq_union_consensus_sv_bedpe
    path hq_union_consensus_upset_intersection_plot

    when:
    params.manta == "on" && params.svaba == "on" && params.delly == "on" && params.igcaller == "on"

    script:
    hq_union_consensus_sv_bedpe = "${tumor_normal_sample_id}.hq.union.consensus.somatic.sv.bedpe"
    hq_union_consensus_upset_intersection_plot = "${tumor_normal_sample_id}.hq.union.consensus.somatic.sv.intersection.plot.pdf"
    """
    Rscript --vanilla ${workflow.projectDir}/bin/sv_union_consensus_polisher.R \
    "${tumor_normal_sample_id}" \
    "${manta_filtered_final_sv_vcf}" \
    "${svaba_filtered_final_sv_vcf}" \
    "${delly_filtered_final_sv_vcf}" \
    "${igcaller_oncogenic_rearrangements_tsv}"
    """
}

// END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \\


// ~~~~~~~~~~~~~~ fragCounter ~~~~~~~~~~~~~~ \\
// START

// fragCounter ~ GC and mappability corrected fragment coverage across a genome for CNV and SV support
process binReadCoverageInNormal_fragcounter {
    publishDir "${params.output_dir}/somatic/fragCounter", mode: 'copy', pattern: '*.{rds,bw}'
    tag "${tumor_normal_sample_id} S=Normal"

    input:
    tuple val(tumor_normal_sample_id), path(normal_bam), path(normal_bam_index) from normal_bams_forFragCounter
    path gc_mappability from gc_mappability_dir
     
    output:
    path normal_fragcounter_cov_rds
    path normal_fragcounter_cov_bw

    when:
    params.fragcounter == "on"

    script:
    normal_id = "${normal_bam.baseName}".replaceFirst(/\..*$/, "")
    normal_results_dir = "results_normal"
    normal_fragcounter_cov_rds = "${normal_id}.fragcounter.cov.rds"
    normal_fragcounter_cov_bw = "${normal_id}.fragcounter.cov.corrected.bw"
    """
    fragCounter_executor.sh \
	"${normal_bam}" \
	"${gc_mappability}" \
	"TRUE" \
	"FALSE" \
	"${normal_results_dir}"

    mv "${normal_results_dir}/cov.rds" "${normal_fragcounter_cov_rds}"
    mv "${normal_results_dir}/cov.corrected.bw" "${normal_fragcounter_cov_bw}"
    """
}

process binReadCoverageInTumor_fragcounter {
    publishDir "${params.output_dir}/somatic/fragCounter", mode: 'copy', pattern: '*.{rds,bw}'
    tag "${tumor_normal_sample_id} S=Tumor"

    input:
    tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index) from tumor_bams_forFragCounter
    path gc_mappability from gc_mappability_dir
     
    output:
    path tumor_fragcounter_cov_rds
    path tumor_fragcounter_cov_bw

    when:
    params.fragcounter == "on"

    script:
    tumor_id = "${tumor_bam.baseName}".replaceFirst(/\..*$/, "")
    tumor_results_dir = "results_tumor"
    tumor_fragcounter_cov_rds = "${tumor_id}.fragcounter.cov.rds"
    tumor_fragcounter_cov_bw = "${tumor_id}.fragcounter.cov.corrected.bw"
    """
    fragCounter_executor.sh \
	"${tumor_bam}" \
	"${gc_mappability}" \
	"TRUE" \
	"FALSE" \
	"${tumor_results_dir}"

    mv "${tumor_results_dir}/cov.rds" "${tumor_fragcounter_cov_rds}"
    mv "${tumor_results_dir}/cov.corrected.bw" "${tumor_fragcounter_cov_bw}"
    """
}

// SNP-Pileup ~ generate SNP read count pileups for CNV calling
process snpPileup_fragcounter {
    publishDir "${params.output_dir}/somatic/fragCounter", mode: 'copy', pattern: '*.{txt.gz}'
    tag "${tumor_normal_sample_id}"

    input:
    tuple val(tumor_normal_sample_id), path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index) from bams_forFragCounterPileup
    path hapmap_vcf from hapmap_ref_snps_vcf
	path hapmap_vcf_index from hapmap_ref_snps_vcf_index

    output:
    tuple val(tumor_normal_sample_id), path(fragcounter_snp_pileup)

    when:
    params.fragcounter == "on"

    script:
    tumor_id = "${tumor_bam.baseName}".replaceFirst(/\..*$/, "")
    normal_id = "${normal_bam.baseName}".replaceFirst(/\..*$/, "")
    tumor_normal_sample_id = "${tumor_id}_vs_${normal_id}"
    fragcounter_snp_pileup = "${tumor_normal_sample_id}.snp_pileup.txt.gz"
    """
    snp-pileup \
    --min-map-quality 0 \
    --min-base-quality 20 \
    --min-read-counts 3 \
    "${hapmap_vcf}" \
    "${tumor_normal_sample_id}.snp_pileup.csv" \
    "${tumor_bam}" \
    "${normal_bam}"

    echo -e "seqname\tstart\tend\talt.count.t\tref.count.t\talt.count.n\tref.count.n" > "${tumor_normal_sample_id}.snp_pileup.txt"
    grep -v 'Chromosome' "${tumor_normal_sample_id}.snp_pileup.csv" \
    | \
    awk -F ',' 'BEGIN {OFS="\t"} {print \$1,\$2,\$2,\$6,\$5,\$10,\$9}' >> "${tumor_normal_sample_id}.snp_pileup.txt"

    gzip "${tumor_normal_sample_id}.snp_pileup.txt"
    """
}

// END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \\


// ~~~~~~~~~~~~~~~ Telomeres ~~~~~~~~~~~~~~~ \\
// START

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
	tuple path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index) from tumor_normal_pair_forTelomereHunter
	path cytoband from cytoband_bed

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
	--bandingFile "${cytoband}" \
	--parallel \
	--plotFileFormat png
	"""
}

// END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \\


