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
log.info "~~~ Run Report File ~~~ 	nextflow_report.${params.run_id}.html"
log.info ''
log.info '################################################'
log.info ''

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
		--help                        [flag]  Prints this message

	Toolbox Switches:
		--telseq                       [str]  Indicates whether or not to use this tool
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
		--controlfreec                 [str]  Indicates whether or not to use this tool
		                                      Available: off, on
		                                      Default: on
		--manta                        [str]  Indicates whether or not to use this tool
		                                      Available: off, on
		                                      Default: on
		--svaba                        [str]  Indicates whether or not to use this tool
		                                      Available: off, on
		                                      Default: on

	################################################

	""".stripIndent()
}

// #################################################### \\
// ~~~~~~~~~~~~~ PARAMETER CONFIGURATION ~~~~~~~~~~~~~~ \\

// Declare the defaults for all pipeline parameters
params.output_dir = "output"
params.run_id = null
params.sample_sheet = null
params.mutect_ref_vcf_concatenated = "yes"
params.vep_ref_cached = "yes"
params.telseq = "on"
params.conpair = "on"
params.varscan = "on"
params.mutect = "on"
params.controlfreec = "on"
params.manta = "on"
params.svaba = "on"
params.help = null

// Print help message if requested
if( params.help ) exit 0, helpMessage()

// Print error messages if required parameters are not set
if( params.run_id == null ) exit 1, "The run command issued does not have the '--run_id' parameter set. Please set the '--run_id' parameter to a unique identifier for the run."

if( params.sample_sheet == null ) exit 1, "The run command issued does not have the '--sample_sheet' parameter set. Please set the '--sample_sheet' parameter to the path of the normal/tumor pair sample sheet CSV."

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
	       reference_genome_fasta_forManta;
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
	       reference_genome_fasta_index_forManta;
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
	       reference_genome_fasta_dict_forManta;
	       reference_genome_fasta_dict_forAnnotation }

Channel
	.fromPath( 'references/hg38/wgs_calling_regions.hg38.bed' )
	.into{ gatk_bundle_wgs_bed_forVarscanSamtoolsMpileup;
	       gatk_bundle_wgs_bed_forMutectCalling;
	       gatk_bundle_wgs_bed_forMutectPileup;
	       gatk_bundle_wgs_bed_forControlFreecSamtoolsMpileup;
	       gatk_bundle_wgs_bed_forManta }

Channel
	.fromList( ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6',
	            'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
	            'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18',
	            'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY',] )
	.into{ chromosome_list_forVarscanSamtoolsMpileup;
	       chromosome_list_forMutectCalling;
	       chromosome_list_forMutectPileup;
	       chromosome_list_forControlFreecSamtoolsMpileup;
	       chromosome_list_forControlFreecMerge }

Channel
	.fromPath( 'references/hg38/sex_identification_loci.chrY.hg38.txt' )
	.set{ sex_identification_loci }

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

if( params.mutect_ref_vcf_concatenated == "yes" ) {
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

if( params.vep_ref_cached == "yes" ) {
	Channel
		.fromPath( 'references/hg38/homo_sapiens_vep_101_GRCh38/', type: 'dir', checkIfExists: true )
		.ifEmpty{ error "The run command issued has the '--vep_ref_cached' parameter set to 'yes', however the directory does not exist. Please set the '--vep_ref_cached' parameter to 'no' and resubmit the run command. For more information, check the README or issue the command 'nextflow run somatic.nf --help'"}
		.set{ vep_ref_dir_preDownloaded }
}

// #################################################### \\
// ~~~~~~~~~~~~~~~~ PIPELINE PROCESSES ~~~~~~~~~~~~~~~~ \\

// Read user provided sample sheet to set the Tumor/Normal sample pairs
Channel
	.fromPath( params.sample_sheet )
	.splitCsv( header:true )
	.map{ row -> def tumor_bam = "input/preprocessedBams/${row.tumor}"
				 def tumor_bam_index = "input/preprocessedBams/${row.tumor}".replaceFirst(/\.bam$/, ".bai")
	             def normal_bam = "input/preprocessedBams/${row.normal}"
	             def normal_bam_index = "input/preprocessedBams/${row.normal}".replaceFirst(/\.bam$/, ".bai")
	             return[ file(tumor_bam), file(tumor_bam_index), file(normal_bam),  file(normal_bam_index) ] }
	.into{ tumor_normal_pair_forAlleleCount;
		   tumor_normal_pair_forConpairPileup;
	       tumor_normal_pair_forVarscanSamtoolsMpileup; 
	       tumor_normal_pair_forVarscanBamReadcount;
	       tumor_normal_pair_forMutectCalling;
	       tumor_normal_pair_forMutectPileup;
	       tumor_normal_pair_forControlFreecSamtoolsMpileup;
	       tumor_normal_pair_forManta }

// Read user provided sample sheet to find Tumor sample BAM files
Channel
	.fromPath( params.sample_sheet )
	.splitCsv( header:true )
	.map{ row -> file("input/preprocessedBams/${row.tumor}") }
	.unique()
	.set{ tumor_bam_forTelomereLengthEstimation }

// Combine reference FASTA index and sex identification loci files into one channel for use in alleleCount process
reference_genome_fasta_index_forAlleleCount.combine( sex_identification_loci )
	.set{ ref_index_and_sex_ident_loci }

// alleleCount ~ determine the sex of each sample to use in downstream analyses
process identifySampleSex_allelecount {
	tag "T=${tumor_id} N=${normal_id}"

	input:
	tuple path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(reference_genome_fasta_index_forAlleleCount), path(sex_identification_loci) from tumor_normal_pair_forAlleleCount.combine(ref_index_and_sex_ident_loci)

	output:
	tuple val(tumor_normal_sample_id), path(sample_sex) into sex_of_sample_forControlFreecCalling

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

// TelSeq ~ estimate telomere length of sample
process tumorSampleTelomereLengthEstimation_telseq {
	publishDir "${params.output_dir}/somatic/telomereLengthEstimations", mode: 'move'
	tag "${tumor_bam.baseName}"

	input:
	path tumor_bam from tumor_bam_forTelomereLengthEstimation

	output:
	path telomere_length_estimation

	when:
	params.telseq == "on"

	script:
	telomere_length_estimation = "${tumor_bam}".replaceFirst(/\..*bam/, ".tumor.telomerelength.txt")
	"""
	telseq "${tumor_bam}" > "${telomere_length_estimation}"
	"""
}

// ~~~~~~~~~~~~~~~~ Conpair ~~~~~~~~~~~~~~~~ \\
// START

// Combine all reference FASTA files into one channel for use in Conpair Pileup process
reference_genome_fasta_forConpairPileup.combine( reference_genome_fasta_index_forConpairPileup )
	.combine( reference_genome_fasta_dict_forConpairPileup )
	.set{ reference_genome_bundle_forConpairPileup }

// Conpair ~ generate GATK pileups the tumor and normal BAMs separately
process bamPileupForConpair_conpair {
	publishDir "${params.output_dir}/somatic/concordanceAndContaminationEstimations", mode: 'symlink'
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

// Conpair ~ concordance and contamination estimator for tumor–normal pileups
process concordanceAndContaminationEstimation_conpair {
	publishDir "${params.output_dir}/somatic/concordanceAndContaminationEstimations", mode: 'move'
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
	publishDir "${params.output_dir}/somatic/varscan", mode: 'symlink'
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
	"${tumor_bam}" "${normal_bam}" \
	| \
	java -jar \${VARSCAN} somatic \
	--mpileup \
	--min-coverage-normal 8 \
	--min-coverage-tumor 6 \
	--min-var-freq 0.08 \
	--min-freq-for-hom 0.75 \
	--normal-purity 1.00 \
	--tumor-purity 1.00 \
	--p-value 0.99 \
	--somatic-p-value 0.05 \
	--strand-filter 1 \
	--output-vcf \
	--output-snp "${tumor_normal_sample_id}.${chromosome}.snv" \
	--output-indel "${tumor_normal_sample_id}.${chromosome}.indel"

	bgzip < "${tumor_normal_sample_id}.${chromosome}.snv.vcf" > "${raw_per_chromosome_snv_vcf}"
	tabix "${raw_per_chromosome_snv_vcf}"

	bgzip < "${tumor_normal_sample_id}.${chromosome}.indel.vcf" > "${raw_per_chromosome_indel_vcf}"
	tabix "${raw_per_chromosome_indel_vcf}"
	"""
}

// BCFtools Concat ~ concatenate all VarScan SNV/indel per chromosome VCFs
process concatenateVarscanPerChromosomeVcfs_bcftools {
	publishDir "${params.output_dir}/somatic/varscan", mode: 'symlink'
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
	publishDir "${params.output_dir}/somatic/varscan", mode: 'symlink'
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

// Combine all needed reference FASTA and input Tumor BAM files into one channel for use in bam-readcount process
reference_genome_fasta_forVarscanBamReadcount.combine( reference_genome_fasta_index_forVarscanBamReadcount )
	.combine( reference_genome_fasta_dict_forVarscanBamReadcount )
	.set{ reference_genome_bundle_forVarscanBamReadcount }

tumor_normal_pair_forVarscanBamReadcount.combine( reference_genome_bundle_forVarscanBamReadcount )
	.set{ bam_and_reference_genome_bundle_forVarscanBamReadcount }

// bam-readcount ~ generate metrics at single nucleotide positions for filtering out false positive calls
process bamReadcountForVarscanFpFilter_bamreadcount {
	publishDir "${params.output_dir}/somatic/varscan", mode: 'symlink'
	tag "${tumor_normal_sample_id}"

	input:
	tuple path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(reference_genome_fasta_forVarscanBamReadcount), path(reference_genome_fasta_index_forVarscanBamReadcount), path(reference_genome_fasta_dict_forVarscanBamReadcount), val(tumor_normal_sample_id), path(high_confidence_snv_vcf), path(high_confidence_snv_vcf_index), path(high_confidence_indel_vcf), path(high_confidence_indel_vcf_index) from bam_and_reference_genome_bundle_forVarscanBamReadcount.combine(high_confidence_vcfs_forVarscanBamReadcount)

	output:
	tuple path(snv_readcount_file), path(indel_readcount_file) into readcount_forVarscanFpFilter

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
	--output "${tumor_normal_sample_id}.somatic.hc.filtered.vcf.gz" \
	"${high_confidence_snv_vcf}" "${high_confidence_indel_vcf}"

	tabix "${tumor_normal_sample_id}.somatic.hc.filtered.vcf.gz"

	bam_readcount_helper.py \
	"${tumor_normal_sample_id}.somatic.hc.filtered.vcf.gz" \
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
	publishDir "${params.output_dir}/somatic/varscan", mode: 'symlink'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(high_confidence_snv_vcf), path(high_confidence_snv_vcf_index), path(high_confidence_indel_vcf), path(high_confidence_indel_vcf_index), path(snv_readcount_file), path(indel_readcount_file) from high_confidence_vcfs_forVarscanFpFilter.combine(readcount_forVarscanFpFilter)
	
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
	--output-file "${fp_filtered_snv_vcf}"

	gunzip -f "${high_confidence_indel_vcf}"

	java -jar \$VARSCAN fpfilter \
	"${unzipped_hc_indel_vcf}" \
	"${indel_readcount_file}" \
	--filtered-file "${tumor_normal_sample_id}.indel.failed.vcf" \
	--output-file "${fp_filtered_indel_vcf}"
	"""
}

// Combine all needed reference FASTA files into one channel for use in VarScan / BCFtools Norm process
reference_genome_fasta_forVarscanBcftoolsNorm.combine( reference_genome_fasta_index_forVarscanBcftoolsNorm )
	.combine( reference_genome_fasta_dict_forVarscanBcftoolsNorm )
	.set{ reference_genome_bundle_forVarscanBcftoolsNorm }

// BCFtools Concat / Norm ~ join the SNV and indel calls, split multiallelic sites into multiple rows then left-align and normalize indels
process concatSplitMultiallelicAndLeftNormalizeVarscanVcf_bcftools {
	publishDir "${params.output_dir}/somatic/varscan", mode: 'symlink'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(fp_filtered_snv_vcf), path(fp_filtered_indel_vcf), path(reference_genome_fasta_forVarscanBcftoolsNorm), path(reference_genome_fasta_index_forVarscanBcftoolsNorm), path(reference_genome_fasta_dict_forVarscanBcftoolsNorm) from filtered_vcfs_forVarscanBcftools.combine(reference_genome_bundle_forVarscanBcftoolsNorm)

	output:
	path final_varscan_vcf into final_varscan_vcf_forAnnotation
	path final_varscan_vcf_index into final_varscan_vcf_index_forAnnotation
	path varscan_multiallelics_stats
	path varscan_realign_normalize_stats

	when:
	params.varscan == "on"

	script:
	final_varscan_vcf = "${tumor_normal_sample_id}.varscan.somatic.vcf.gz"
	final_varscan_vcf_index = "${final_varscan_vcf}.tbi"
	varscan_multiallelics_stats = "${tumor_normal_sample_id}.varscan.multiallelicsstats.txt"
	varscan_realign_normalize_stats = "${tumor_normal_sample_id}.varscan.realignnormalizestats.txt"
	"""
	bcftools concat \
	--threads ${task.cpus} \
	--allow-overlaps \
	--output-type v \
	"${fp_filtered_snv_vcf}" "${fp_filtered_indel_vcf}" \
	| \
	bgzip --stdout \
	| \
	bcftools norm \
	--threads ${task.cpus} \
	--multiallelics -both \
	--output-type z \
	- 2>"${varscan_multiallelics_stats}" \
	| \
	bcftools norm \
	--threads ${task.cpus} \
	--fasta-ref "${reference_genome_fasta_forVarscanBcftoolsNorm}" \
	--output-type z \
	- 2>"${varscan_realign_normalize_stats}" \
	--output "${final_varscan_vcf}"

	tabix "${final_varscan_vcf}"
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
if( params.mutect_ref_vcf_concatenated == "yes" ) {
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
process mutect2SnvAndIndelCalling_gatk {
	publishDir "${params.output_dir}/somatic/mutect", mode: 'symlink'
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
	publishDir "${params.output_dir}/somatic", mode: 'symlink'
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
	publishDir "${params.output_dir}/somatic/mutect", mode: 'symlink'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(raw_per_chromosome_mutect_stats_file) from raw_per_chromosome_mutect_stats_forMutectStatsMerge.groupTuple()

	output:
	path merged_mutect_stats_file into merged_mutect_stats_file_forMutectFilter

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
	publishDir "${params.output_dir}/somatic/mutect", mode: 'symlink'
	tag "C=${chromosome} T=${tumor_id} N=${normal_id}"

	input:
	tuple path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(gatk_bundle_wgs_bed_forMutectPileup), path(exac_common_sites_ref_vcf), path(exac_common_sites_ref_vcf_index) from tumor_normal_pair_forMutectPileup.combine(bed_and_resources_vcfs_forMutectPileup)
	each chromosome from chromosome_list_forMutectPileup

	output:
	tuple val(tumor_id), path(per_chromosome_tumor_pileup) into per_chromosome_tumor_pileups_forMutectPileupGather
	tuple val(normal_id), path(per_chromosome_normal_pileup) into per_chromosome_normal_pileups_forMutectPileupGather

	when:
	params.mutect == "on"

	script:
	tumor_id = "${tumor_bam.baseName}".replaceFirst(/\..*$/, "")
	normal_id = "${normal_bam.baseName}".replaceFirst(/\..*$/, "")
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
	publishDir "${params.output_dir}/somatic/mutect", mode: 'symlink'
	tag "${tumor_id}"

	input:
	tuple val(tumor_id), path(per_chromosome_tumor_pileup), path(reference_genome_fasta_dict) from per_chromosome_tumor_pileups_forMutectPileupGather.groupTuple().combine(reference_genome_fasta_dict_forMutectPileupGatherTumor)

	output:
	tuple val(tumor_id), path(tumor_pileup) into tumor_pileups_forMutectContamination

	when:
	params.mutect == "on"

	script:
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
	publishDir "${params.output_dir}/somatic/mutect", mode: 'symlink'
	tag "${normal_id}"

	input:
	tuple val(normal_id), path(per_chromosome_normal_pileup), path(reference_genome_fasta_dict) from per_chromosome_normal_pileups_forMutectPileupGather.groupTuple().combine(reference_genome_fasta_dict_forMutectPileupGatherNormal)

	output:
	tuple val(normal_id), path(normal_pileup) into normal_pileups_forMutectContamination

	when:
	params.mutect == "on"

	script:
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
	publishDir "${params.output_dir}/somatic/mutect", mode: 'symlink'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_id), path(tumor_pileup), val(normal_id), path(normal_pileup) from tumor_pileups_forMutectContamination.combine(normal_pileups_forMutectContamination)

	output:
	path contamination_file into contamination_file_forMutectFilter

	when:
	params.mutect == "on"

	script:
	tumor_normal_sample_id = "${tumor_id}_vs_${normal_id}"
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

merged_raw_vcfs_forMutectFilter.combine( merged_mutect_stats_file_forMutectFilter )
	.combine( contamination_file_forMutectFilter )
	.set{ input_vcf_stats_and_contamination_forMutectFilter }

// GATK FilterMutectCalls ~ filter somatic SNVs and indels called by Mutect2
process mutect2VariantFiltration_gatk {
	publishDir "${params.output_dir}/somatic/mutect", mode: 'symlink'
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
	publishDir "${params.output_dir}/somatic/mutect", mode: 'symlink'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(filtered_vcf), path(filtered_vcf_index), path(reference_genome_fasta_forMutectBcftools), path(reference_genome_fasta_index_forMutectBcftools), path(reference_genome_fasta_dict_forMutectBcftools) from filtered_vcf_forMutectBcftools.combine(reference_genome_bundle_forMutectBcftools)

	output:
	path final_mutect_vcf into final_mutect_vcf_forAnnotation
	path final_mutect_vcf_index into final_mutect_vcf_index_forAnnotation
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
	publishDir "${params.output_dir}/somatic/controlFreec", mode: 'symlink'
	tag "C=${chromosome} T=${tumor_id} N=${normal_id}"

	input:
	tuple path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(reference_genome_fasta_forControlFreecSamtoolsMpileup), path(reference_genome_fasta_index_forControlFreecSamtoolsMpileup), path(reference_genome_fasta_dict_forControlFreecSamtoolsMpileup), path(gatk_bundle_wgs_bed_forControlFreecSamtoolsMpileup) from tumor_normal_pair_forControlFreecSamtoolsMpileup.combine(reference_genome_bundle_and_bed_forControlFreecSamtoolsMpileup)
	each chromosome from chromosome_list_forControlFreecSamtoolsMpileup

	output:
	tuple val(tumor_id), path(tumor_pileup_per_chromosome) into per_chromosome_tumor_pileups_forControlFreecMerge
	tuple val(normal_id), path(normal_pileup_per_chromosome) into per_chromosome_normal_pileups_forControlFreecMerge

	when:
	params.controlfreec == "on"

	script:
	tumor_id = "${tumor_bam.baseName}".replaceFirst(/\..*$/, "")
	normal_id = "${normal_bam.baseName}".replaceFirst(/\..*$/, "")
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
	publishDir "${params.output_dir}/somatic/controlFreec", mode: 'symlink'
	tag "T=${tumor_id} N=${normal_id}"

	input:
	tuple val(tumor_id), path(tumor_pileup_per_chromosome) from per_chromosome_tumor_pileups_forControlFreecMerge.groupTuple()
	tuple val(normal_id), path(normal_pileup_per_chromosome) from per_chromosome_normal_pileups_forControlFreecMerge.groupTuple()
	val chromosome_list from chromosome_list_forControlFreecMerge.collect()

	output:
	tuple val(tumor_normal_sample_id), path(tumor_pileup), path(normal_pileup) into tumor_normal_pileups_forControlFreecCalling

	when:
	params.controlfreec == "on"

	script:
	tumor_normal_sample_id = "${tumor_id}_vs_${normal_id}"
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
tumor_normal_pileups_forControlFreecCalling.groupTuple( sex_of_sample_forControlFreecCalling )
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
	publishDir "${params.output_dir}/somatic/controlFreec", mode: 'symlink'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(tumor_pileup), path(normal_pileup), path(sex_of_sample_forControlFreecCalling), path(reference_genome_fasta_forControlFreecCalling), path(reference_genome_fasta_index_forControlFreecCalling), path(reference_genome_fasta_dict_forControlFreecCalling), path(autosome_sex_chromosome_fasta_dir), path(autosome_sex_chromosome_sizes), path(common_dbsnp_ref_vcf), path(common_dbsnp_ref_vcf_index), path(mappability_track_zip) from tumor_normal_pileups_and_sex_ident.combine(reference_data_bundle_forControlFreec)

	//tuple val(tumor_normal_sample_id), path(tumor_pileup), path(normal_pileup), path(reference_genome_fasta_forControlFreecCalling), path(reference_genome_fasta_index_forControlFreecCalling), path(reference_genome_fasta_dict_forControlFreecCalling), path(autosome_sex_chromosome_fasta_dir), path(autosome_sex_chromosome_sizes), path(common_dbsnp_ref_vcf), path(common_dbsnp_ref_vcf_index), path(mappability_track_zip) from tumor_normal_pileups_forControlFreecCalling.combine(reference_data_bundle_forControlFreec)

	output:
	tuple val(tumor_normal_sample_id), path(cnv_profile_raw), path(cnv_ratio_file), path(baf_file) into cnv_calling_files_forControlFreecPostProcessing
	path config_file
	path subclones_file

	when:
	params.controlfreec == "on"

	script:
	config_file = "${tumor_normal_sample_id}.controlfreec.config.txt"
	cnv_profile_raw = "${tumor_normal_sample_id}.controlfreec.raw.cnv"
	cnv_ratio_file = "${tumor_normal_sample_id}.controlfreec.ratio.txt"
	subclones_file = "${tumor_normal_sample_id}.controlfreec.subclones.txt"
	baf_file = "${tumor_normal_sample_id}.controlfreec.baf.txt"
	"""
	unzip -q "${mappability_track_zip}"
	gunzip -q "${autosome_sex_chromosome_fasta_dir}"/*.fa.gz
	sex=\$(cut -d ' ' -f 2 "${sex_of_sample_forControlFreecCalling}")

	touch "${config_file}"
	echo "[general]" >> "${config_file}"
	echo "chrFiles = \${PWD}/${autosome_sex_chromosome_fasta_dir}" >> "${config_file}"
	echo "chrLenFile = \${PWD}/${autosome_sex_chromosome_sizes}" >> "${config_file}"
	echo "gemMappabilityFile = \${PWD}/out100m2_hg38.gem" >> "${config_file}"
	echo "minimalSubclonePresence = 20" >> "${config_file}"
	echo "maxThreads = ${task.cpus}" >> "${config_file}"
	echo "ploidy = 2" "${config_file}"
	echo "sex = \${sex}" >> "${config_file}"
	echo "window = 50000" >> "${config_file}"
	echo "" >> "${config_file}"

	echo "[sample]" >> "${config_file}"
	echo "mateFile = ${tumor_pileup}" >> "${config_file}"
	echo "inputFormat = pileup" >> "${config_file}"
	echo "mateOrientation = FR" >> "${config_file}"
	echo "" >> "${config_file}"

	echo "[control]" >> "${config_file}"
	echo "mateFile = ${normal_pileup}" >> "${config_file}"
	echo "inputFormat = pileup" >> "${config_file}"
	echo "mateOrientation = FR" >> "${config_file}"
	echo "" >> "${config_file}"

	echo "[BAF]" >> "${config_file}"
	echo "SNPfile = \${PWD}/${common_dbsnp_ref_vcf}" >> "${config_file}"
	echo "" >> "${config_file}"

	freec -conf "${config_file}"

	mv "${tumor_pileup}_CNVs" "${cnv_profile_raw}"
	mv "${tumor_pileup}_ratio.txt" "${cnv_ratio_file}"
	mv "${tumor_pileup}_subclones.txt" "${subclones_file}"
	mv "${tumor_pileup}_BAF.txt" "${baf_file}"
	"""
}

// Control-FREEC ~ post-processing of CNV predictions for significance, visualization, and format compatibility
process cnvPredictionPostProcessing_controlfreec {
	publishDir "${params.output_dir}/somatic/controlFreec", mode: 'symlink'
	tag "${tumor_normal_sample_id}"

	input:
	tuple val(tumor_normal_sample_id), path(cnv_profile_raw), path(cnv_ratio_file), path(baf_file) from cnv_calling_files_forControlFreecPostProcessing

	output:
	path cnv_profile_final
	path cnv_ratio_bed_file

	when:
	params.controlfreec == "on"

	script:
	cnv_profile_final = "${tumor_normal_sample_id}.controlfreec.cnv.txt"
	cnv_ratio_bed_file = "${tumor_normal_sample_id}.controlfreec.ratio.bed"
	"""
	cat \${CONTROLFREEC_DIR}/scripts/assess_significance.R | R --slave --args "${cnv_profile_raw}" "${cnv_ratio_file}"
	mv "${cnv_profile_raw}.p.value.txt" "${cnv_profile_final}"

	cat \${CONTROLFREEC_DIR}/scripts/makeGraph.R | R --slave --args 2 "${cnv_ratio_file}" "${baf_file}"


	perl \${CONTROLFREEC_DIR}/scripts/freec2bed.pl -f "${cnv_ratio_file}" > "${cnv_ratio_bed_file}"
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
	publishDir "${params.output_dir}/somatic/manta", mode: 'symlink'
	tag "T=${tumor_id} N=${normal_id}"

	input:
	tuple path(tumor_bam), path(tumor_bam_index), path(normal_bam), path(normal_bam_index), path(reference_genome_fasta_forManta), path(reference_genome_fasta_index_forManta), path(reference_genome_fasta_dict_forManta), path(gatk_bundle_wgs_bed_forManta) from tumor_normal_pair_forManta.combine(reference_genome_bundle_and_bed_forManta)

	output:
	path final_manta_somatic_sv_vcf into final_manta_vcf_forAnnotation
	path final_manta_somatic_sv_vcf_index into final_manta_vcf_index_forAnnotation
	tuple path(candidate_sv_vcf), path(candidate_sv_vcf_index)
	tuple path(candidate_indels_vcf), path(candidate_indels_vcf_index)
	tuple path(diploid_sv_vcf), path(diploid_sv_vcf_index)

	when:
	params.manta == "on"

	script:
	tumor_id = "${tumor_bam.baseName}".replaceFirst(/\..*$/, "")
	normal_id = "${normal_bam.baseName}".replaceFirst(/\..*$/, "")
	tumor_normal_sample_id = "${tumor_id}_vs_${normal_id}"
	zipped_wgs_bed = "${gatk_bundle_wgs_bed_forManta}.gz"
	candidate_sv_vcf = "${tumor_normal_sample_id}.manta.candidateSV.vcf.gz"
	candidate_sv_vcf_index = "${candidate_sv_vcf}.tbi"
	candidate_indels_vcf = "${tumor_normal_sample_id}.manta.candidateSmallIndels.vcf.gz"
	candidate_indels_vcf_index = "${candidate_indels_vcf}.tbi"
	diploid_sv_vcf = "${tumor_normal_sample_id}.manta.diploid.vcf.gz"
	diploid_sv_vcf_index = "${diploid_sv_vcf}.tbi"
	final_manta_somatic_sv_vcf = "${tumor_normal_sample_id}.manta.somatic.vcf.gz"
	final_manta_somatic_sv_vcf_index = "${final_manta_somatic_sv_vcf}.tbi"
	"""
	bgzip < "${gatk_bundle_wgs_bed_forManta}" > "${zipped_wgs_bed}"
	tabix "${zipped_wgs_bed}"

	python \${MANTA_DIR}/bin/configManta.py \
	--tumorBam "${tumor_bam}" \
	--normalBam "${normal_bam}" \
	--referenceFasta "${reference_genome_fasta_forManta}" \
	--callRegions "${zipped_wgs_bed}" \
	--runDir manta

	python manta/runWorkflow.py \
	--mode local \
	--jobs ${task.cpus} \
	--memGb ${task.memory.toGiga()}

	mv manta/results/variants/candidateSV.vcf.gz "${candidate_sv_vcf}"
	mv manta/results/variants/candidateSV.vcf.gz.tbi "${candidate_sv_vcf_index}"
	mv manta/results/variants/candidateSmallIndels.vcf.gz "${candidate_indels_vcf}"
	mv manta/results/variants/candidateSmallIndels.vcf.gz.tbi "${candidate_indels_vcf_index}"

	zcat manta/results/variants/diploidSV.vcf.gz \
	| \
	grep -E "^#|PASS" \
	| \
	bgzip > "${diploid_sv_vcf}"
	tabix "${diploid_sv_vcf}"

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



/*


// ~~~~~~~~~~~~~~~~~ SvABA ~~~~~~~~~~~~~~~~~ \\
// START

// Combine
BWArefgenomefilesALL
HomosapiansdbSNPvcf
wgsBed?


// SvABA ~ detecting structural variants using genome-wide local assembly
process svAndIndelCalling_svaba {
	publishDir "${params.output_dir}/somatic/svaba", mode: 'symlink'
	tag "T=${tumor_id} N=${normal_id}"

	input:


	output:


	when:
	params.svaba == "on"

	script:

	"""
	svaba run \
	-t "${TUMBAM}" \
	-n "${NORMBAM}" \
	--reference-genome "${BWAREFGENOMEFILES}" \
	## POTENTIALLY --region "${WGSBED}" \ 
	--id-string "${TUMNORMID}" \
	--dbsnp-vcf "${HOMOSAPDBSNPVCF}" \
	--threads "${task.cpus}" \
	--verbose 1 \
	--g-zip
	"""
}

// END
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \\


*/




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

// Create a channel of each unannotated somatic VCF and index produced in the pipeline to put through the VEP annotation process
final_varscan_vcf_forAnnotation.ifEmpty{ 'skipped' }
	.mix( final_mutect_vcf_forAnnotation.ifEmpty{ 'skipped' } )
	.mix( final_manta_vcf_forAnnotation.ifEmpty{ 'skipped' } )
	.filter{ it != 'skipped' }
	.set{ final_somatic_vcfs_forAnnotation }

final_varscan_vcf_index_forAnnotation.ifEmpty{ 'skipped' }
	.mix( final_mutect_vcf_index_forAnnotation.ifEmpty{ 'skipped' } )
	.mix( final_manta_vcf_index_forAnnotation.ifEmpty{ 'skipped' } )
	.filter{ it != 'skipped' }
	.set{ final_somatic_vcfs_indicies_forAnnotation }

// VEP ~ annotate the final somatic VCFs using databases including Ensembl, GENCODE, RefSeq, PolyPhen, SIFT, dbSNP, COSMIC, etc.
process annotateSomaticVcf_vep {
	publishDir "${params.output_dir}/somatic/vepAnnotatedVcfs", mode: 'copy'
	tag "${vcf_id}"

	input:
	tuple path(cached_ref_dir_vep), path(reference_genome_fasta_forAnnotation), path(reference_genome_fasta_index_forAnnotation), path(reference_genome_fasta_dict_forAnnotation) from vep_ref_dir.combine(reference_genome_bundle_forAnnotation)
	path final_somatic_vcfs_indicies_forAnnotation
	each path(final_somatic_vcf) from final_somatic_vcfs_forAnnotation

	output:
	path final_annotated_somatic_vcfs
	path annotation_summary

	script:
	vcf_id = "${final_somatic_vcf}".replaceFirst(/\.somatic\.vcf\.gz/, "")
	final_annotated_somatic_vcfs = "${vcf_id}.annotated.somatic.vcf.gz"
	annotation_summary = "${final_somatic_vcf}".replaceFirst(/\.vcf\.gz/, ".vep.summary.html")
	"""
	vep \
	--offline \
	--cache \
	--dir "${cached_ref_dir_vep}" \
	--assembly GRCh38 \
	--fasta "${reference_genome_fasta_forAnnotation}" \
	--input_file "${final_somatic_vcf}" \
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
