// Myeloma Genome Project 1000
// Comprehensive pipeline for analysis of matched T/N Multiple Myeloma WGS data
// https://github.com/pblaney/mgp1000

// This portion of the pipeline is used for germline variant analysis of normal WGS samples.
// It is designed to be run with BAMs that were genereated via the Preprocessing step of
// this pipeline.

import java.text.SimpleDateFormat;
def workflowTimestamp = "${workflow.start.format('MM-dd-yyyy HH:mm')}"

log.info ''
log.info '##### Myeloma Genome Project 1000 Pipeline #####'
log.info '################################################'
log.info '~~~~~~~~~~~ GERMLINE VARIANT ANALYSIS ~~~~~~~~~~'
log.info '################################################'
log.info ''
log.info "~~~ Launch Time ~~~		${workflowTimestamp}"
log.info ''
log.info "~~~ Output Directory ~~~ 	${workflow.projectDir}/output/germline"
log.info ''
log.info '################################################'
log.info ''

def helpMessage() {
	log.info"""

	Usage Example:

		nextflow run germline.nf -bg -resume --sample_sheet samplepairlist.csv --cohort_name batch1 --singularity_module singularity/3.1 --email someperson@gmail.com --vep_ref_cached no --ref_vcf_concatenated no -profile germline 

	Mandatory Arguments:
    	--sample_sheet                 [str]  CSV file containing the list of samples where the first column designates the file name of the
    	                                      normal sample, the second column for the file name of the match tumor sample, and the third
    	                                      column has the sex of the sample, example of the format for this file is in the testSamples
    	                                      directory
    	--cohort_name                  [str]  A user defined collective name of the group of samples being run through this step of the
    	                                      pipeline, this will be used as the name of the final output genotyped VCF
		--email                        [str]  Email address to send workflow completion/stoppage notification
		-profile                       [str]  Configuration profile to use, each profile described in nextflow.config file
		                                      Currently available: preprocessing

	Main Options:
		-bg                           [flag]  Runs the pipeline processes in the background, this option should be included if deploying
		                                      pipeline with real data set so processes will not be cut if user disconnects from deployment
		                                      environment
		-resume                       [flag]  Successfully completed tasks are cached so that if the pipeline stops prematurely the
		                                      previously completed tasks are skipped while maintaining their output
		--singularity_module           [str]  Indicates the name of the Singularity software module to be loaded for use in the pipeline,
		                                      this option is not needed if Singularity is natively installed on the deployment environment
		--vep_ref_cached               [str]  Indicates whether or not the VEP reference files used for annotation have been downloaded/cached
		                                      locally, this will be done in a process of the pipeline if it has not, this does not need to be
		                                      done for every separate run after the first
		--ref_vcf_concatenated         [str]  Indicates whether or not the 1000 Genomes Project reference VCF used for ADMIXTURE analysis has
		                                      been concatenated, this will be done in a process of the pipeline if it has not, this does not
		                                      need to be done for every separate run after the first 
		--help                        [flag]  Prints this message

	################################################

	""".stripIndent()
}

// #################################################### \\
// ~~~~~~~~~~~~~ PARAMETER CONFIGURATION ~~~~~~~~~~~~~~ \\

// Declare the defaults for all pipeline parameters
params.output_dir = "output"
params.sample_sheet = null
params.cohort_name = null
params.vep_ref_cached = "yes"
params.ref_vcf_concatenated = "yes"
params.help = null

// Print help message if requested
if( params.help ) exit 0, helpMessage()

// Set channels for reference files
Channel
	.fromPath( 'references/hg38/Homo_sapiens_assembly38.fasta' )
	.into{ reference_genome_fasta_forSplitIntervals;
		   reference_genome_fasta_forHaplotypeCaller;
		   reference_genome_fasta_forCombineGvcfs;
	       reference_genome_fasta_forJointGenotyping;
	       reference_genome_fasta_forIndelVariantRecalibration;
	       reference_genome_fasta_forSnpVariantRecalibration;
	       reference_genome_fasta_forSplitAndNorm;
	       reference_genome_fasta_forAnnotation }

Channel
	.fromPath( 'references/hg38/Homo_sapiens_assembly38.fasta.fai' )
	.into{ reference_genome_fasta_index_forSplitIntervals;
	       reference_genome_fasta_index_forHaplotypeCaller;
	       reference_genome_fasta_index_forCombineGvcfs;
		   reference_genome_fasta_index_forJointGenotyping;
		   reference_genome_fasta_index_forIndelVariantRecalibration;
		   reference_genome_fasta_index_forSnpVariantRecalibration;
		   reference_genome_fasta_index_forSplitAndNorm;
		   reference_genome_fasta_index_forAnnotation }

Channel
	.fromPath( 'references/hg38/Homo_sapiens_assembly38.dict' )
	.into{ reference_genome_fasta_dict_forSplitIntervals;
	       reference_genome_fasta_dict_forHaplotypeCaller;
	       reference_genome_fasta_dict_forCombineGvcfs;
	       reference_genome_fasta_dict_forJointGenotyping;
	       reference_genome_fasta_dict_forIndelVariantRecalibration;
	       reference_genome_fasta_dict_forSnpVariantRecalibration;
	       reference_genome_fasta_dict_forSplitAndNorm;
	       reference_genome_fasta_dict_forAnnotation }

Channel
	.fromPath( 'references/hg38/wgs_calling_regions.hg38.interval_list' )
	.set{ gatk_bundle_wgs_interval_list }

Channel
	.fromPath( 'references/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz' )
	.set{ gatk_bundle_mills_1000G }

Channel
	.fromPath( 'references/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi' )
	.set{ gatk_bundle_mills_1000G_index }

Channel
	.fromPath( 'references/hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz' )
	.into{ gatk_bundle_dbsnp138_forJointGenotyping;
	       gatk_bundle_dbsnp138_forIndelVariantRecalibration;
	       gatk_bundle_dbsnp138_forSnpVariantRecalibration }

Channel
	.fromPath( 'references/hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi' )
	.into{ gatk_bundle_dbsnp138_index_forJointGenotyping;
	       gatk_bundle_dbsnp138_index_forIndelVariantRecalibration;
	       gatk_bundle_dbsnp138_index_forSnpVariantRecalibration }

Channel
	.fromPath( 'references/hg38/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz' )
	.set{ gatk_bundle_axiom }

Channel
	.fromPath( 'references/hg38/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi' )
	.set{ gatk_bundle_axiom_index }

Channel
	.fromPath( 'references/hg38/Hapmap_3.3.hg38.vcf.gz' )
	.set{ gatk_bundle_hapmap }

Channel
	.fromPath( 'references/hg38/Hapmap_3.3.hg38.vcf.gz.tbi' )
	.set{ gatk_bundle_hapmap_index }

Channel
	.fromPath( 'references/hg38/1000G_Omni2.5.hg38.vcf.gz' )
	.set{ gatk_bundle_1000G_omni }

Channel
	.fromPath( 'references/hg38/1000G_Omni2.5.hg38.vcf.gz.tbi' )
	.set{ gatk_bundle_1000G_omni_index }

Channel
	.fromPath( 'references/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz' )
	.set{ gatk_bundle_1000G_snps }

Channel
	.fromPath( 'references/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi' )
	.set{ gatk_bundle_1000G_snps_index }

if( params.vep_ref_cached == "yes" ) {
	Channel
		.fromPath( 'references/hg38/homo_sapiens_vep_101_GRCh38/', type: dir, checkIfExists: true )
		.ifEmpty{ error "The run command issued has the '--vep_ref_cached' parameter set to 'yes', however the directory does not exist. Please set the '--vep_ref_cached' parameter to 'no' and resubmit the run command. For more information, check the README or issue the command 'nextflow run germline.nf --help'"}
		.set{ vep_ref_dir_preDownloaded }
}

Channel
	.fromPath( 'references/hg38/ALL.chr[0-9].shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz' )
	.toSortedList()
	.set{ reference_vcf_1000G_chromosomes1_9 }

Channel
	.fromPath( 'references/hg38/ALL.chr[0-9][0-9].shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz' )
	.toSortedList()
	.set{ reference_vcf_1000G_chromosomes10_22 }

Channel
	.fromPath( 'references/hg38/ALL.chrX.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz' )
	.set{ reference_vcf_1000G_chromosomeX }

Channel
	.fromPath( 'references/hg38/ALL.chr*.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz.tbi' )
	.set{ reference_vcf_1000G_per_chromosome_index }

Channel
	.fromPath( 'references/hg38/ref_vcf_samples_to_keep.txt' )
	.set{ reference_vcf_1000G_samples_to_keep }

Channel
	.fromPath( 'references/hg38/ref_vcf_chr_name_conversion_map.txt' )
	.set{ reference_vcf_1000G_chr_name_conversion_map }

Channel
	.fromPath( 'references/hg38/ref_vcf_samples_known_ancestry.txt' )
	.set{ reference_vcf_1000G_known_ancestry }

if( params.ref_vcf_concatenated == "yes" ) {
	Channel
		.fromPath( 'references/hg38/ALL.wgs.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz', checkIfExists: true )
		.ifEmpty{ error "The run command issued has the '--ref_vcf_concatenated' parameter set to 'yes', however the file does not exist. Please set the '--ref_vcf_concatenated' parameter to 'no' and resubmit the run command. For more information, check the README or issue the command 'nextflow run germline.nf --help'"}
		.set{ reference_vcf_1000G_preBuilt }

	Channel
		.fromPath( 'references/hg38/ALL.wgs.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz.tbi', checkIfExists: true )
		.ifEmpty{ error "The '--ref_vcf_concatenated' parameter set to 'yes', however the index file does not exist for the reference VCF. Please set the '--ref_vcf_concatenated' parameter to 'no' and resubmit the run command. Alternatively, use Tabix to index the reference VCF."}
		.set{ reference_vcf_1000G_index_preBuilt }
}


// #################################################### \\
// ~~~~~~~~~~~~~~~~ PIPELINE PROCESSES ~~~~~~~~~~~~~~~~ \\

// Read user provided sample sheet to find Normal sample BAM files
Channel
	.fromPath( params.sample_sheet )
	.ifEmpty{ error "No sample sheet provided, an example is given in the testSamples directory" }
	.splitCsv( header:true )
	.map{ row -> file("input/preprocessedBams/${row.normal}") }
	.unique()
	.into{ input_preprocessed_bams_forTelomereLengthEstimation;
	       input_preprocessed_bams_forHaplotypeCaller }

// TelSeq ~ estimate telomere length of sample
process telomereLengthEstimation_telseq {
	publishDir "${params.output_dir}/germline/telomereLengthEstimations", mode: 'move'
	tag "${bam_preprocessed.baseName}"

	input:
	path bam_preprocessed from input_preprocessed_bams_forTelomereLengthEstimation

	output:
	path telomere_length_estimation

	script:
	telomere_length_estimation = "${bam_preprocessed}".replaceFirst(/\..*bam/, ".telomerelength.csv")
	"""
	telseq "${bam_preprocessed}" > "${telomere_length_estimation}"
	"""
}

// Combine all needed reference FASTA files into one channel for use in SplitIntervals process
reference_genome_fasta_forSplitIntervals.combine( reference_genome_fasta_index_forSplitIntervals )
	.combine( reference_genome_fasta_dict_forSplitIntervals )
	.set{ reference_genome_bundle_forSplitIntervals }

// GATK SplitIntervals ~ divide interval list into files containing equal number of reads for better pipeline performance
process splitIntervalList_gatk {
	tag "Spliting interval list into 20 files"
	
	input:
	tuple path(reference_genome_fasta_forSplitIntervals), path(reference_genome_fasta_index_forSplitIntervals), path(reference_genome_fasta_dict_forSplitIntervals) from reference_genome_bundle_forSplitIntervals
	path gatk_bundle_wgs_interval_list

	output:
	path "splitIntervals/*-split.interval_list" into split_intervals mode flatten

	script:
	"""
	gatk SplitIntervals \
	--java-options "-Xmx${task.memory.toGiga()}G -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1 -XX:+AggressiveOpts" \
	--verbosity ERROR \
	--tmp-dir . \
	--reference "${reference_genome_fasta_forSplitIntervals}" \
	--intervals "${gatk_bundle_wgs_interval_list}" \
	--subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION \
	--extension -split.interval_list \
	--scatter-count 20 \
	--output splitIntervals
	"""
}

// Combine all needed reference FASTA files into one channel for use in GATK HaplotypeCaller process
reference_genome_fasta_forHaplotypeCaller.combine( reference_genome_fasta_index_forHaplotypeCaller )
	.combine( reference_genome_fasta_dict_forHaplotypeCaller )
	.set{ reference_genome_bundle_forHaplotypeCaller }

// Combine the input BAM file channel with the reference FASTA channel
input_preprocessed_bams_forHaplotypeCaller.combine( reference_genome_bundle_forHaplotypeCaller )
	.set{ input_bams_and_reference_fasta_forHaplotypeCaller }

// GATK HaplotypeCaller ~ call germline SNPs and indels via local re-assembly
process haplotypeCaller_gatk {
	tag "${sample_id}.${interval_id}"

	input:
	tuple path(bam_preprocessed), path(reference_genome_fasta_forHaplotypeCaller), path(reference_genome_fasta_index_forHaplotypeCaller), path(reference_genome_fasta_dict_forHaplotypeCaller), path(interval) from input_bams_and_reference_fasta_forHaplotypeCaller.combine(split_intervals)

	output:
	tuple val(sample_id), path(gvcf_per_interval_raw), path(gvcf_per_interval_raw_index) into raw_gvcfs

	script:
	sample_id = "${bam_preprocessed}".replaceFirst(/\.final\..*bam/, "")
	interval_id = "${interval}".replaceFirst(/-split\.interval_list/, "_int")
	gvcf_per_interval_raw = "${sample_id}.${interval_id}.g.vcf"
	gvcf_per_interval_raw_index = "${gvcf_per_interval_raw}.idx"
	"""
	gatk HaplotypeCaller \
	--java-options "-Xmx${task.memory.toGiga()}G -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1 -XX:+AggressiveOpts" \
	--verbosity ERROR \
	--tmp-dir . \
	--max-alternate-alleles 3 \
	--standard-min-confidence-threshold-for-calling 50 \
	--emit-ref-confidence GVCF \
	--reference "${reference_genome_fasta_forHaplotypeCaller}" \
	--intervals "${interval}" \
	--input "${bam_preprocessed}" \
	--output "${gvcf_per_interval_raw}"
	"""
}

// GATK SortVcfs ~ merge all GVCF files for each sample and sort them
process mergeAndSortGvcfs_gatk {
	publishDir "${params.output_dir}/germline/mergedAndSortedGvcfs", mode: 'symlink'
	tag "${sample_id}"
	
	input:
	tuple val(sample_id), path(gvcf_per_interval_raw), path(gvcf_per_interval_raw_index) from raw_gvcfs.groupTuple()

	output:
	path gvcf_merged_raw into merged_raw_gcvfs
	path gvcf_merged_raw_index

	script:
	gvcf_merged_raw = "${sample_id}.g.vcf"
	gvcf_merged_raw_index = "${gvcf_merged_raw}.idx"
	"""
	gatk SortVcf \
	--java-options "-Xmx${task.memory.toGiga()}G -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1 -XX:+AggressiveOpts" \
	--VERBOSITY ERROR \
	--TMP_DIR . \
	--MAX_RECORDS_IN_RAM 4000000 \
	${gvcf_per_interval_raw.collect { "--INPUT $it " }.join()} \
	--OUTPUT "${gvcf_merged_raw}"
	"""
}

// Combine all needed reference FASTA files into one channel for use in GATK CombineGVCFs process
reference_genome_fasta_forCombineGvcfs.combine( reference_genome_fasta_index_forCombineGvcfs)
	.combine( reference_genome_fasta_dict_forCombineGvcfs )
	.set{ reference_genome_bundle_forCombineGvcfs }

// GATK CombineGVFCs ~ combine per-sample GVCF files into a multi-sample GVCF for joint-genotyping
process combineAllGvcfs_gatk {
	publishDir "${params.output_dir}/germline/combinedGvcf", mode: 'symlink'
	tag "${params.cohort_name} GVCF"

	input:
	path gvcf_merged_raw from merged_raw_gcvfs.toList()
	tuple path(reference_genome_fasta_forCombineGvcfs), path(reference_genome_fasta_index_forCombineGvcfs), path(reference_genome_fasta_dict_forCombineGvcfs) from reference_genome_bundle_forCombineGvcfs

	output:
	tuple path(gvcf_cohort_combined), path(gvcf_cohort_combined_index) into combined_cohort_gvcf

	script:
	gvcf_cohort_combined_unzipped = "${params.cohort_name}.g.vcf"
	gvcf_cohort_combined = "${params.cohort_name}.g.vcf.gz"
	gvcf_cohort_combined_index = "${gvcf_cohort_combined}.tbi"
	"""
	gatk CombineGVCFs \
	--java-options "-Xmx${task.memory.toGiga()}G -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1 -XX:+AggressiveOpts" \
	--verbosity ERROR \
	--tmp-dir . \
	--reference "${reference_genome_fasta_forCombineGvcfs}" \
	${gvcf_merged_raw.collect {" --variant $it" }.join()} \
	--output "${gvcf_cohort_combined_unzipped}"

	bgzip < "${gvcf_cohort_combined_unzipped}" > "${gvcf_cohort_combined}"
	tabix "${gvcf_cohort_combined}"
	"""
}

// Combine all needed GATK bundle files and reference FASTA files into one channel for use in GATK GenotypeGVCFs process
gatk_bundle_dbsnp138_forJointGenotyping.combine( gatk_bundle_dbsnp138_index_forJointGenotyping )
	.set{ gatk_reference_bundle_forJointGenotyping }

reference_genome_fasta_forJointGenotyping.combine( reference_genome_fasta_index_forJointGenotyping )
	.combine( reference_genome_fasta_dict_forJointGenotyping )
	.set{ reference_genome_bundle_forJointGenotyping }

// GATK GenotypeGVCFs ~ perform joint genotyping
process jointGenotyping_gatk {
	publishDir "${params.output_dir}/germline/genotypedVcf", mode: 'symlink'

	input:
	tuple path(gvcf_cohort_combined), path(gvcf_cohort_combined_index) from combined_cohort_gvcf
	tuple path(reference_genome_fasta_forJointGenotyping), path(reference_genome_fasta_index_forJointGenotyping), path(reference_genome_fasta_dict_forJointGenotyping) from reference_genome_bundle_forJointGenotyping
	tuple path(gatk_bundle_dbsnp138), path(gatk_bundle_dbsnp138_index) from gatk_reference_bundle_forJointGenotyping

	output:
	tuple path(vcf_joint_genotyped), path(vcf_joint_genotyped_index) into joint_genotyped_vcfs

	script:
	vcf_joint_genotyped = "${params.cohort_name}.vcf.gz"
	vcf_joint_genotyped_index = "${vcf_joint_genotyped}.tbi"
	"""
	gatk GenotypeGVCFs \
	--java-options "-Xmx${task.memory.toGiga()}G -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1 -XX:+AggressiveOpts" \
	--verbosity ERROR \
	--tmp-dir . \
	--dbsnp "${gatk_bundle_dbsnp138}" \
	--reference "${reference_genome_fasta_forJointGenotyping}" \
	--standard-min-confidence-threshold-for-calling 50 \
	--annotation-group StandardAnnotation \
	--annotation-group AS_StandardAnnotation \
	--variant "${gvcf_cohort_combined}" \
	--output "${vcf_joint_genotyped}"
	"""
}

// GATK VariantFiltration ~ hard filter test for excess heterozygosity then remove non-site level genotype information
// ExcessHet is a phred-scaled p-value. We want a cutoff of anything more extreme
// than a z-score of -4.5 which is a p-value of 3.4e-06, which phred-scaled is 54.69
process excessHeterozygosityHardFilter_gatk {
	publishDir "${params.output_dir}/germline/excessHetHardFilter", mode: 'symlink'
	tag "Excess Heterozygosity Hard Filter and Non-Site Genotype Removal"

	input:
	tuple path(vcf_joint_genotyped), path(vcf_joint_genotyped_index) from joint_genotyped_vcfs

	output:
	tuple path(vcf_hard_filtered), path(vcf_hard_filtered_index) into hard_filtered_vcfs_forIndelVariantRecalibration, hard_filtered_vcfs_forSnvVariantRecalibration, hard_filtered_vcfs_forApplyVqsr

	script:
	vcf_hard_filtered = "${vcf_joint_genotyped}".replaceFirst(/\.vcf\.gz/, ".filtered.vcf.gz")
	vcf_hard_filtered_index = "${vcf_hard_filtered}.tbi"
	"""
	gatk VariantFiltration \
	--java-options "-Xmx${task.memory.toGiga()}G -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1 -XX:+AggressiveOpts" \
	--verbosity ERROR \
	--tmp-dir . \
	--filter-name ExcessHet \
	--filter-expression "ExcessHet > 54.69" \
	--variant "${vcf_joint_genotyped}" \
	--output "${vcf_hard_filtered}"
	"""
}

// Combine all needed GATK bundle files and reference FASTA files into one channel for use in GATK Indel VariantRecalibratior process
gatk_bundle_mills_1000G.combine( gatk_bundle_mills_1000G_index )
	.combine( gatk_bundle_axiom )
	.combine( gatk_bundle_axiom_index )
	.combine( gatk_bundle_dbsnp138_forIndelVariantRecalibration )
	.combine( gatk_bundle_dbsnp138_index_forIndelVariantRecalibration)
	.set{ gatk_reference_bundle_forIndelVariantRecalibration}

reference_genome_fasta_forIndelVariantRecalibration.combine( reference_genome_fasta_index_forIndelVariantRecalibration )
	.combine( reference_genome_fasta_dict_forIndelVariantRecalibration )
	.set{ reference_genome_bundle_forIndelVariantRecalibration }

// GATK VariantRecalibrator (Indels) ~ build recalibration model to score indel variant quality for filtering
process indelVariantRecalibration_gatk {
	publishDir "${params.output_dir}/germline/indelVariantRecal", mode: 'symlink'
	tag "Indel Variant Recalibration"

	input:
	tuple path(vcf_hard_filtered), path(vcf_hard_filtered_index) from hard_filtered_vcfs_forIndelVariantRecalibration
	tuple path(gatk_bundle_mills_1000G), path(gatk_bundle_mills_1000G_index), path(gatk_bundle_axiom), path(gatk_bundle_axiom_index), path(gatk_bundle_dbsnp138_forIndelVariantRecalibration), path(gatk_bundle_dbsnp138_index_forIndelVariantRecalibration) from gatk_reference_bundle_forIndelVariantRecalibration
	tuple path(reference_genome_fasta_forIndelVariantRecalibration), path(reference_genome_fasta_index_forIndelVariantRecalibration), path(reference_genome_fasta_dict_forIndelVariantRecalibration) from reference_genome_bundle_forIndelVariantRecalibration

	output:
	tuple path(indel_vqsr_table), path(indel_vqsr_table_index), path(indel_vqsr_tranches) into indel_vqsr_files

	script:
	indel_vqsr_table = "${params.cohort_name}.indel.recaldata.table"
	indel_vqsr_table_index = "${indel_vqsr_table}.idx"
	indel_vqsr_tranches = "${params.cohort_name}.indel.tranches"
	"""
	gatk VariantRecalibrator \
	--java-options "-Xmx${task.memory.toGiga()}G -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1 -XX:+AggressiveOpts" \
	--verbosity ERROR \
	--tmp-dir . \
	--mode INDEL \
	-an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
	--max-gaussians 4 \
	--trust-all-polymorphic \
	--reference "${reference_genome_fasta_forIndelVariantRecalibration}" \
	--resource:mills,known=false,training=true,truth=true,prior=12 "${gatk_bundle_mills_1000G}" \
	--resource:axiomPoly,known=false,training=true,truth=false,prior=10 "${gatk_bundle_axiom}" \
	--resource:dbsnp,known=true,training=false,truth=false,prior=2 "${gatk_bundle_dbsnp138_forIndelVariantRecalibration}" \
	--variant "${vcf_hard_filtered}" \
	--output "${indel_vqsr_table}" \
	--tranches-file "${indel_vqsr_tranches}"
	"""
}

// Combine all needed GATK bundle files and reference FASTA files into one channel for use in GATK SNV VariantRecalibratior process
gatk_bundle_hapmap.combine( gatk_bundle_hapmap_index )
	.combine( gatk_bundle_1000G_omni )
	.combine( gatk_bundle_1000G_omni_index )
	.combine( gatk_bundle_1000G_snps )
	.combine( gatk_bundle_1000G_snps_index)
	.combine( gatk_bundle_dbsnp138_forSnpVariantRecalibration )
	.combine( gatk_bundle_dbsnp138_index_forSnpVariantRecalibration )
	.set{ gatk_reference_bundle_forSnvVariantRecalibration }

reference_genome_fasta_forSnpVariantRecalibration.combine( reference_genome_fasta_index_forSnpVariantRecalibration )
	.combine( reference_genome_fasta_dict_forSnpVariantRecalibration )
	.set{ reference_genome_bundle_forSnpVariantRecalibration }

// GATK VariantRecalibrator (SNPs) ~ build recalibration model to score SNP variant quality for filtering
process snpVariantRecalibration_gatk {
	publishDir "${params.output_dir}/germline/snpVariantRecal", mode: 'symlink'
	tag "SNP Variant Recalibration"

	input:
	tuple path(vcf_hard_filtered), path(vcf_hard_filtered_index) from hard_filtered_vcfs_forSnvVariantRecalibration
	tuple path(gatk_bundle_hapmap), path(gatk_bundle_hapmap_index), path(gatk_bundle_1000G_omni), path(gatk_bundle_1000G_omni_index), path(gatk_bundle_1000G_snps), path(gatk_bundle_1000G_snps_index), path(gatk_bundle_dbsnp138_forSnpVariantRecalibration), path(gatk_bundle_dbsnp138_index_forSnpVariantRecalibration) from gatk_reference_bundle_forSnvVariantRecalibration
	tuple path(reference_genome_fasta_forSnpVariantRecalibration), path(reference_genome_fasta_index_forSnpVariantRecalibration), path(reference_genome_fasta_dict_forSnpVariantRecalibration) from reference_genome_bundle_forSnpVariantRecalibration

	output:
	tuple path(snp_vqsr_table), path(snp_vqsr_table_index), path(snp_vqsr_tranches) into snp_vqsr_files

	script:
	snp_vqsr_table = "${params.cohort_name}.snp.recaldata.table"
	snp_vqsr_table_index = "${snp_vqsr_table}.idx"
	snp_vqsr_tranches = "${params.cohort_name}.snp.tranches"	
	"""
	gatk VariantRecalibrator \
	--java-options "-Xmx${task.memory.toGiga()}G -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1 -XX:+AggressiveOpts" \
	--verbosity ERROR \
	--tmp-dir . \
	--mode SNP \
	-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP \
	--max-gaussians 6 \
	--trust-all-polymorphic \
	--reference "${reference_genome_fasta_forSnpVariantRecalibration}" \
	--resource:hapmap,known=false,training=true,truth=true,prior=15 "${gatk_bundle_hapmap}" \
	--resource:omni,known=false,training=true,truth=true,prior=12 "${gatk_bundle_1000G_omni}" \
	--resource:1000G,known=false,training=true,truth=false,prior=10 "${gatk_bundle_1000G_snps}" \
	--resource:dbsnp,known=true,training=false,truth=false,prior=7 "${gatk_bundle_dbsnp138_forSnpVariantRecalibration}" \
	--variant "${vcf_hard_filtered}" \
	--output "${snp_vqsr_table}" \
	--tranches-file "${snp_vqsr_tranches}"
	"""
}

// GATK ApplyVQSR ~ apply variant quality score recalibration for Indels and SNPs
process applyIndelAndSnpVqsr_gatk {
	publishDir "${params.output_dir}/germline/vqsrVcfs", mode: 'symlink'
	tag "Applying Variant Quality Score Recalibration to Indels and SNPs"

	input:
	tuple path(vcf_hard_filtered), path(vcf_hard_filtered_index) from hard_filtered_vcfs_forApplyVqsr
	tuple path(indel_vqsr_table), path(indel_vqsr_table_index), path(indel_vqsr_tranches) from indel_vqsr_files
	tuple path(snp_vqsr_table), path(snp_vqsr_table_index), path(snp_vqsr_tranches) from snp_vqsr_files

	output:
	tuple path(final_vqsr_germline_vcf), path(final_vqsr_germline_vcf_index) into vqsr_germline_vcfs

	script:
	intermediate_vqsr_germline_vcf = "${vcf_hard_filtered}".replaceFirst(/\.filtered\.vcf\.gz/, ".intermediate.vqsr.vcf.gz")
	final_vqsr_germline_vcf = "${vcf_hard_filtered}".replaceFirst(/\.filtered\.vcf\.gz/, ".final.vqsr.vcf.gz")
	final_vqsr_germline_vcf_index = "${final_vqsr_germline_vcf}.tbi"
	"""
	gatk ApplyVQSR \
	--java-options "-Xmx${task.memory.toGiga()}G -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1 -XX:+AggressiveOpts" \
	--verbosity ERROR \
	--tmp-dir . \
	--mode INDEL \
	--truth-sensitivity-filter-level 99.0 \
	--exclude-filtered \
	--variant "${vcf_hard_filtered}" \
	--recal-file "${indel_vqsr_table}" \
	--tranches-file "${indel_vqsr_tranches}" \
	--create-output-variant-index true \
	--output "${intermediate_vqsr_germline_vcf}"

	gatk ApplyVQSR \
	--java-options "-Xmx${task.memory.toGiga()}G -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1 -XX:+AggressiveOpts" \
	--verbosity ERROR \
	--tmp-dir . \
	--mode SNP \
	--truth-sensitivity-filter-level 99.5 \
	--exclude-filtered \
	--variant "${intermediate_vqsr_germline_vcf}" \
	--recal-file "${snp_vqsr_table}" \
	--tranches-file "${snp_vqsr_tranches}" \
	--create-output-variant-index true \
	--output "${final_vqsr_germline_vcf}"
	"""
}

//Combine all needed reference FASTA files into one channel for use in BCFtools Norm process
reference_genome_fasta_forSplitAndNorm.combine( reference_genome_fasta_index_forSplitAndNorm )
	.combine( reference_genome_fasta_dict_forSplitAndNorm )
	.set{ reference_genome_bundle_forSplitAndNorm }

// BCFtools Norm ~ Split multiallelic sites into multiple rows then left-align and normalize indels
process splitMultiallelicAndLeftNormalizeVcf_bcftools {
	publishDir "${params.output_dir}/germline/finalGermlineVcfs", mode: 'symlink'
	tag "Splitting multiallelic sites and left-normalizing indels"

	input:
	tuple path(final_vqsr_germline_vcf), path(final_vqsr_germline_vcf_index) from vqsr_germline_vcfs
	tuple path(reference_genome_fasta_forSplitAndNorm), path(reference_genome_fasta_index_forSplitAndNorm), path(reference_genome_fasta_dict_forSplitAndNorm) from reference_genome_bundle_forSplitAndNorm

	output:
	tuple path(final_germline_vcf), path(final_germline_vcf_index) into final_germline_vcf_forAnnotation, final_germline_vcf_forAdmixture
	path multiallelics_stats
	path realign_normalize_stats

	script:
	final_germline_vcf = "${final_vqsr_germline_vcf}".replaceFirst(/\.final\.vqsr\.vcf\.gz/, ".germline.vcf.gz")
	final_germline_vcf_index = "${final_germline_vcf}.tbi"
	multiallelics_stats = "${params.cohort_name}.multiallelicsstats.txt"
	realign_normalize_stats = "${params.cohort_name}.realignnormalizestats.txt"
	"""
	zcat "${final_vqsr_germline_vcf}" \
	| \
	bcftools norm \
	--threads ${task.cpus} \
	--multiallelics -both \
	--output-type z \
	- 2>"${multiallelics_stats}" \
	| \
	bcftools norm \
	--threads ${task.cpus} \
	--fasta-ref "${reference_genome_fasta_forSplitAndNorm}" \
	--output-type z \
	- 2>"${realign_normalize_stats}" \
	> "${final_germline_vcf}"

	tabix "${final_germline_vcf}"
	"""
}

// VEP ~ Download the reference files used for VEP annotation, if needed
process downloadVepAnnotationReferences_vep {
	publishDir "references/hg38", mode: 'copy'
	tag "Downloading references files for VEP annotation"

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

//Combine all needed reference FASTA files into one channel for use in VEP annotation process
reference_genome_fasta_forAnnotation.combine( reference_genome_fasta_index_forAnnotation )
	.combine( reference_genome_fasta_dict_forAnnotation )
	.set{ reference_genome_bundle_forAnnotation }

// VEP ~ Annotate the final germline VCF using databases including Ensembl, GENCODE, RefSeq, PolyPhen, SIFT, dbSNP, COSMIC, etc.
process annotateGermlineVcf_vep {
	publishDir "${params.output_dir}/germline/vepAnnotatedVcf", mode: 'symlink'
	tag "Annotating ${params.cohort} germline VCF"

	input:
	tuple path(final_germline_vcf), path(final_germline_vcf_index) from final_germline_vcf_forAnnotation
	path cached_ref_dir_vep from vep_ref_dir
	tuple path(reference_genome_fasta_forAnnotation), path(reference_genome_fasta_index_forAnnotation), path(reference_genome_fasta_dict_forAnnotation) from reference_genome_bundle_forAnnotation

	output:
	path final_annotated_germline_vcf
	path annotation_summary

	script:
	final_annotated_germline_vcf = "${final_germline_vcf}".replaceFirst(/\.germline\.vcf\.gz/, "annotated.germline.vcf.gz")
	annotation_summary = "${params.cohort_name}.vep.summary.html"
	"""
	vep \
	--offline \
	--cache \
	--dir "${cached_ref_dir_vep}" \
	--assembly GRCh38 \
	--fasta "${reference_genome_fasta_forAnnotation}" \
	--input_file "${final_germline_vcf}" \
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
	--output_file "${final_annotated_germline_vcf}" \
	--compress_output bgzip \
	--vcf
	"""
}

// Merge all reference VCFs per chromosome into one sorted channel for concatenation
reference_vcf_1000G_chromosomes1_9
	.merge( reference_vcf_1000G_chromosomes10_22 )
	.merge( reference_vcf_1000G_chromosomeX )
	.set{ reference_vcf_1000G_per_chromosome }

// BCFtools Concat/Annotate/View ~ Prepare the 1000 Genomes Project reference VCFs for use in ADMIXTURE process, if needed
process referenceVcfPrep_bcftools {
	publishDir "references/hg38", mode: 'copy'
	tag "Concatenating reference genome for ADMIXTURE projection analysis"

	input:
	path per_chromosome_ref_vcf from reference_vcf_1000G_per_chromosome
	path per_chromosome_ref_vcf_index from reference_vcf_1000G_per_chromosome_index.toList()
	path ref_vcf_samples_to_keep_list from reference_vcf_1000G_samples_to_keep
	path chr_name_conversion_map from reference_vcf_1000G_chr_name_conversion_map

	output:
	tuple path(whole_genome_ref_vcf), path(whole_genome_ref_vcf_index) into reference_vcf_1000G_fromProcess

	when:
	params.ref_vcf_concatenated == "no"

	script:
	whole_genome_ref_vcf = "ALL.wgs.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
	whole_genome_ref_vcf_index = "ALL.wgs.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz.tbi"
	"""
	bcftools concat \
	--threads ${task.cpus} \
	--output-type z \
	${per_chromosome_ref_vcf.collect { "$it "}.join()} \
	| \
	bcftools annotate \
	--threads ${task.cpus} \
	--output-type z \
	--rename-chrs "${chr_name_conversion_map}" \
	- \
	| \
	bcftools view \
	--threads ${task.cpus} \
	--output-type z \
	--samples-file "${ref_vcf_samples_to_keep_list}" \
	- > "${whole_genome_ref_vcf}"

	tabix "${whole_genome_ref_vcf}"
	"""
}

// Depending on whether the 1000 Genomes Project whole-genome reference VCF was pre-built, set the input
// channel for the BCFtools merge VCF process
if( params.ref_vcf_concatenated == "yes" ) {
	reference_vcf_1000G = reference_vcf_1000G_preBuilt.combine( reference_vcf_1000G_index_preBuilt )
}
else {
	reference_vcf_1000G = reference_vcf_1000G_fromProcess
}

// BCFtools Annotate/Merge ~ Preprocessing of VCF to remove all FORMAT fields except genotype needed for merging with reference VCF
// and merge cohort VCF with 1000G ref VCF for supervised projection analysis, output no new multiallelic records
process mergeCohortAndReferenceVcf_bcftools {
	publishDir "${params.output_dir}/germline/cohortAndRefMergedVcf", mode: 'symlink'
	tag "Merging cohort VCF with reference VCF"

	input:
	tuple path(vcf_germline_final), path(vcf_germline_final_index) from final_germline_vcf_forAdmixture
	tuple path(whole_genome_ref_vcf), path(whole_genome_ref_vcf_index) from reference_vcf_1000G

	output:
	path sample_ref_merged_vcf into merged_unfiltered_vcf

	script:
	sample_ref_merged_vcf = "${params.cohort_name}.refmerged.vcf.gz"
	vcf_germline_final_gt_only = "${vcf_germline_final}".replaceFirst(/\.germline\.vcf\.gz/, ".germline.gto.vcf.gz")
	"""
	bcftools annotate \
	--threads ${task.cpus} \
	--remove FORMAT \
	--output-type z \
	--output "${vcf_germline_final_gt_only}" \
	"${vcf_germline_final}"

	tabix "${vcf_germline_final_gt_only}"

	bcftools merge \
	--threads ${task.cpus} \
	--merge none \
	--missing-to-ref \
	--output-type z \
	--output "${sample_ref_merged_vcf}" \
	"${vcf_germline_final_gt_only}" "${whole_genome_ref_vcf}"
	"""
}

// VCFtools ~ Hard filter the merged VCF to only contain biallelic, non-singleton SNP sites that are a minimum of 2kb apart from each other
process hardFilterCohortReferenceMergedVcf_vcftools {
	publishDir "${params.output_dir}/germline/hardFilteredMergedVcfPlinkFiles", mode: 'symlink'
	tag "Hard filtering cohort and reference merged VCF"

	input:
	path sample_ref_merged_vcf from merged_unfiltered_vcf

	output:
	tuple val(hard_filtered_plink_file_prefix), path(hard_filtered_plink_ped_file), path(hard_filtered_plink_map_file) into hard_filtered_refmerged_plink_files

	script:
	hard_filtered_plink_file_prefix = "${params.cohort_name}.hardfiltered.refmerged"
	hard_filtered_plink_ped_file = "${hard_filtered_plink_file_prefix}.ped"
	hard_filtered_plink_map_file = "${hard_filtered_plink_file_prefix}.map"
	"""
	vcftools \
	--gzvcf "${sample_ref_merged_vcf}" \
	--thin 2000 \
	--min-alleles 2 \
	--max-alleles 2 \
	--non-ref-ac 2 \
	--plink \
	--temp . \
	--out "${hard_filtered_plink_file_prefix}"
	"""
}

// PLINK ~ Generate ADMIXTURE ready PLINK files keeping only sites with a minor allele freq > 0.05, no missing genotype, then
// prune the markers for linkage disequilibrium (remove SNPs that have an R-squared value of greater than 0.5 with any
// other SNP within a 50-SNP sliding window, the window is advanced by 10-SNPs each time)
process filterPlinkFilesForAdmixture_plink {
	publishDir "${params.output_dir}/germline/mafGenotypeAndLinkeageDiseqFilteredPlinkFiles", mode: 'symlink'
	tag "Filtering PLINK files for MAF < 0.05, no missing genotypes, and pruned for linkage disequilibrium"

	input:
	tuple val(hard_filtered_plink_file_prefix), path(hard_filtered_plink_ped_file), path(hard_filtered_plink_map_file) from hard_filtered_refmerged_plink_files

	output:
	tuple val(pruned_filtered_plink_file_prefix), path(pruned_filtered_plink_bed_file), path(pruned_filtered_plink_bim_file), path(pruned_filtered_plink_fam_file) into pruned_filtered_refmerged_plink_files

	script:
	maf_gt_filtered_plink_file_prefix = "${params.cohort_name}.maf.gt.filtered.refmerged"
	pruned_filtered_plink_file_prefix = "${params.cohort_name}.pruned.maf.gt.filtered.refmerged"
	pruned_filtered_plink_bed_file = "${pruned_filtered_plink_file_prefix}.bed"
	pruned_filtered_plink_bim_file = "${pruned_filtered_plink_file_prefix}.bim"
	pruned_filtered_plink_fam_file = "${pruned_filtered_plink_file_prefix}.fam"
	"""
	plink \
	--threads ${task.cpus} \
	--file "${hard_filtered_plink_file_prefix}" \
	--out "${maf_gt_filtered_plink_file_prefix}" \
	--make-bed \
	--maf 0.05 \
	--geno 0.1 \
	--indep-pairwise 50 10 0.5

	plink \
	--threads ${task.cpus} \
	--bfile "${maf_gt_filtered_plink_file_prefix}" \
	--extract "${maf_gt_filtered_plink_file_prefix}.prune.in" \
	--make-bed \
	--out "${pruned_filtered_plink_file_prefix}"
	"""
}

// ADMIXTURE ~ estimation of sample ancestry using autosomal SNP genotype data in a supervised and haploid aware fashion 
process ancestryEstimation_admixture {
	publishDir "${params.output_dir}/germline/admixutreAncestryEstimation", mode: 'symlink'
	tag "Calculating ancestry estimations"

	input:
	tuple val(pruned_filtered_plink_file_prefix), path(pruned_filtered_plink_bed_file), path(pruned_filtered_plink_bim_file), path(pruned_filtered_plink_fam_file) from pruned_filtered_refmerged_plink_files
	path known_ancestry_file_ref_vcf from reference_vcf_1000G_known_ancestry

	output:
	path admixture_ancestry_fractions
	path admixture_allele_frequencies
	path admixture_standard_error

	script:
	admixture_ancestry_fractions = "${pruned_filtered_plink_file_prefix}.13.Q"
	admixture_allele_frequencies = "${pruned_filtered_plink_file_prefix}.13.P"
	admixture_standard_error = "${pruned_filtered_plink_file_prefix}.13.Q_se"
	"""
	cohort_pop_file_creator.sh \
	"${pruned_filtered_plink_file_prefix}.fam" \
	"${known_ancestry_file_ref_vcf}" \
	"${pruned_filtered_plink_file_prefix}.pop"

	admixture \
	-j${task.cpus}\
	-B50 \
	--supervised \
	--haploid="male:23" \
	"${pruned_filtered_plink_file_prefix}.bed" \
	13
	"""
}
