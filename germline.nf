// Myeloma Genome Project 1000
// Comprehensive pipeline for analysis of matched T/N Multiple Myeloma WGS data
// https://github.com/pblaney/mgp1000

// This portion of the pipeline is used for germline variant analysis of normal WGS samples.
// It can be run on any preprocessed BAM files, even if not run through the Preprocessing step of
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

		nextflow run germline.nf -bg -resume --singularity_module singularity/3.1 --email someperson@gmail.com -profile germline 

	Mandatory Arguments:
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
		--help                        [flag]  Prints this message

	################################################

	""".stripIndent()
}

// #################################################### \\
// ~~~~~~~~~~~~~ PARAMETER CONFIGURATION ~~~~~~~~~~~~~~ \\

// Declare the defaults for all pipeline parameters
params.output_dir = "output"
params.cohort_name = null
params.sample_sheet = "samplesheet.csv"
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
	       reference_genome_fasta_forSnpVariantRecalibration }

Channel
	.fromPath( 'references/hg38/Homo_sapiens_assembly38.fasta.fai' )
	.into{ reference_genome_fasta_index_forSplitIntervals;
	       reference_genome_fasta_index_forHaplotypeCaller;
	       reference_genome_fasta_index_forCombineGvcfs;
		   reference_genome_fasta_index_forJointGenotyping;
		   reference_genome_fasta_index_forIndelVariantRecalibration;
		   reference_genome_fasta_index_forSnpVariantRecalibration }

Channel
	.fromPath( 'references/hg38/Homo_sapiens_assembly38.dict' )
	.into{ reference_genome_fasta_dict_forSplitIntervals;
	       reference_genome_fasta_dict_forHaplotypeCaller;
	       reference_genome_fasta_dict_forCombineGvcfs;
	       reference_genome_fasta_dict_forJointGenotyping;
	       reference_genome_fasta_dict_forIndelVariantRecalibration;
	       reference_genome_fasta_dict_forSnpVariantRecalibration }

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


// #################################################### \\
// ~~~~~~~~~~~~~~~~ PIPELINE PROCESSES ~~~~~~~~~~~~~~~~ \\

// Read user provided sample sheet to find Normal sample BAM files
Channel
	.fromPath( params.samplesheet )
	.ifEmpty{ error "No sample sheet provided, an example is given in the testSamples directory" }
	.splitCsv( header:true )
	.map{ row -> file("input/preprocessedBams/${row.normal}") }
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
	gvcf_cohort_combined = "${params.cohort_name}.g.vcf"
	gvcf_cohort_combined_index = "${gvcf_cohort_combined}.idx"
	"""
	gatk CombineGVCFs \
	--java-options "-Xmx${task.memory.toGiga()}G -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1 -XX:+AggressiveOpts" \
	--verbosity ERROR \
	--tmp-dir . \
	--reference "${reference_genome_fasta_forCombineGvcfs}" \
	${gvcf_merged_raw.collect {" --variant $it" }.join()} \
	--output "${gvcf_cohort_combined}"
	"""
}

// Combine all needed GATK bundle files and reference FASTA files into one channel for use in GATK GenotypeGVCFs process
gatk_bundle_dbsnp138_forJointGenotyping.combine( gatk_bundle_dbsnp138_index_forJointGenotyping )
	.set{ gatk_reference_bundle_forJointGenotyping }

reference_genome_fasta_forJointGenotyping.combine( reference_genome_fasta_index_forJointGenotyping )
	.combine( reference_genome_fasta_dict_forJointGenotyping )
	.set{ reference_genome_bundle_forJointGenotyping }

// GATK GenotypeGVCFs ~ perform joint genotyping using the GenomicsDB workspace
process jointGenotyping_gatk {
	publishDir "${params.output_dir}/germline/genotypedVcf", mode: 'symlink'

	input:
	tuple path(gvcf_cohort_combined), path(gvcf_cohort_combined_index) from combined_cohort_gvcf
	tuple path(reference_genome_fasta_forJointGenotyping), path(reference_genome_fasta_index_forJointGenotyping), path(reference_genome_fasta_dict_forJointGenotyping) from reference_genome_bundle_forJointGenotyping
	tuple path(gatk_bundle_dbsnp138), path(gatk_bundle_dbsnp138_index) from gatk_reference_bundle_forJointGenotyping

	output:
	tuple path(vcf_joint_genotyped), path(vcf_joint_genotyped_index) into joint_genotyped_vcfs

	script:
	vcf_joint_genotyped = "${params.cohort_name}.vcf"
	vcf_joint_genotyped_index = "${vcf_joint_genotyped}.idx"
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

// GATK VariantFiltration/MakeSitesOnlyVcf ~ hard filter test for excess heterozygosity then remove non-site level genotype information
// ExcessHet is a phred-scaled p-value. We want a cutoff of anything more extreme
// than a z-score of -4.5 which is a p-value of 3.4e-06, which phred-scaled is 54.69
process excessHeterozygosityHardFilter_gatk {
	publishDir "${params.output_dir}/germline/excessHetHardFilter", mode: 'symlink'
	tag "Excess Heterozygosity Hard Filter and Non-Site Genotype Removal"

	input:
	tuple path(vcf_joint_genotyped), path(vcf_joint_genotyped_index) from joint_genotyped_vcfs

	output:
	tuple path(vcf_filtered_sites_only), path(vcf_filtered_sites_only_index) into filtered_sites_only_vcfs_forIndelVariantRecalibration, filtered_sites_only_vcfs_forSnvVariantRecalibration, filtered_sites_only_vcfs_forApplyVqsr

	script:
	intermediate_vcf_genotyped_filtered = "${vcf_joint_genotyped}".replaceFirst(/\.vcf/, ".markedforfilter.vcf")
	vcf_filtered_sites_only = "${vcf_joint_genotyped}".replaceFirst(/\.vcf/, ".filtered.vcf")
	vcf_filtered_sites_only_index = "${vcf_filtered_sites_only}.idx"
	"""
	gatk VariantFiltration \
	--java-options "-Xmx${task.memory.toGiga()}G -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1 -XX:+AggressiveOpts" \
	--verbosity ERROR \
	--tmp-dir . \
	--filter-name ExcessHet \
	--filter-expression "ExcessHet > 54.69" \
	--variant "${vcf_joint_genotyped}" \
	--output "${intermediate_vcf_genotyped_filtered}"

	gatk MakeSitesOnlyVcf \
	--java-options "-Xmx${task.memory.toGiga()}G -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1 -XX:+AggressiveOpts" \
	--VERBOSITY ERROR \
	--TMP_DIR . \
	--INPUT "${intermediate_vcf_genotyped_filtered}" \
	--OUTPUT "${vcf_filtered_sites_only}"
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
	tuple path(vcf_filtered_sites_only), path(vcf_filtered_sites_only_index) from filtered_sites_only_vcfs_forIndelVariantRecalibration
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
	--variant "${vcf_filtered_sites_only}" \
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
	tuple path(vcf_filtered_sites_only), path(vcf_filtered_sites_only_index) from filtered_sites_only_vcfs_forSnvVariantRecalibration
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
	--variant "${vcf_filtered_sites_only}" \
	--output "${snp_vqsr_table}" \
	--tranches-file "${snp_vqsr_tranches}"
	"""
}

// GATK ApplyVQSR ~ apply variant quality score recalibration for Indels and SNPs
process applyIndelAndSnpVqsr_gatk {
	publishDir "${params.output_dir}/germline/finalGermlineVcfs", mode: 'symlink'
	tag "Applying Variant Quality Score Recalibration to Indels and SNPs"

	input:
	tuple path(vcf_filtered_sites_only), path(vcf_filtered_sites_only_index) from filtered_sites_only_vcfs_forApplyVqsr
	tuple path(indel_vqsr_table), path(indel_vqsr_table_index), path(indel_vqsr_tranches) from indel_vqsr_files
	tuple path(snp_vqsr_table), path(snp_vqsr_table_index), path(snp_vqsr_tranches) from snp_vqsr_files

	output:
	tuple path(vcf_germline_final), path(vcf_germline_final_index) into final_germline_vcfs

	script:
	intermediate_vcf_germline_final = "${vcf_filtered_sites_only}".replaceFirst(/\.filtered\.vcf/, ".intermediate.vcf")
	vcf_germline_final = "${vcf_filtered_sites_only}".replaceFirst(/\.filtered\.vcf/, ".final.vcf")
	vcf_germline_final_index = "${vcf_germline_final}.idx"
	"""
	gatk ApplyVQSR \
	--java-options "-Xmx${task.memory.toGiga()}G -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1 -XX:+AggressiveOpts" \
	--verbosity ERROR \
	--tmp-dir . \
	--mode INDEL \
	--truth-sensitivity-filter-level 99.0 \
	--exclude-filtered \
	--variant "${vcf_filtered_sites_only}" \
	--recal-file "${indel_vqsr_table}" \
	--tranches-file "${indel_vqsr_tranches}" \
	--create-output-variant-index true \
	--output "${intermediate_vcf_germline_final}"

	gatk ApplyVQSR \
	--java-options "-Xmx${task.memory.toGiga()}G -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1 -XX:+AggressiveOpts" \
	--verbosity ERROR \
	--tmp-dir . \
	--mode SNP \
	--truth-sensitivity-filter-level 99.5 \
	--exclude-filtered \
	--variant "${intermediate_vcf_germline_final}" \
	--recal-file "${snp_vqsr_table}" \
	--tranches-file "${snp_vqsr_tranches}" \
	--create-output-variant-index true \
	--output "${vcf_germline_final}"
	"""
}
