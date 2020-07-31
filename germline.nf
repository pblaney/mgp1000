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
	.into{ reference_genome_fasta_forHaplotypeCaller;
	       reference_genome_fasta_forJointGenotyping;
	       reference_genome_fasta_forIndelVariantRecalibration;
	       reference_genome_fasta_forSnpVariantRecalibration }

Channel
	.fromPath( 'references/hg38/Homo_sapiens_assembly38.fasta.fai' )
	.into{ reference_genome_fasta_index_forHaplotypeCaller;
		   reference_genome_fasta_index_forJointGenotyping;
		   reference_genome_fasta_index_forIndelVariantRecalibration;
		   reference_genome_fasta_index_forSnpVariantRecalibration }

Channel
	.fromPath( 'references/hg38/Homo_sapiens_assembly38.dict' )
	.into{ reference_genome_fasta_dict_forHaplotypeCaller;
	       reference_genome_fasta_dict_forJointGenotyping;
	       reference_genome_fasta_dict_forIndelVariantRecalibration;
	       reference_genome_fasta_dict_forSnpVariantRecalibration }

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

// Combine all needed reference FASTA files into one channel for use in GATK HaplotypeCaller process
reference_genome_fasta_forHaplotypeCaller.combine( reference_genome_fasta_index_forHaplotypeCaller )
	.combine( reference_genome_fasta_dict_forHaplotypeCaller )
	.set{ reference_genome_bundle_forHaplotypeCaller }

// GATK HaplotypeCaller ~ call germline SNPs and indels via local re-assembly
process haplotypeCaller_gatk {
	publishDir "${params.output_dir}/germline/haplotypeCaller", mode: 'symlink'
	tag "${bam_preprocessed.baseName}"

	input:
	tuple path(bam_preprocessed), path(reference_genome_fasta_forHaplotypeCaller), path(reference_genome_fasta_index_forHaplotypeCaller), path(reference_genome_fasta_dict_forHaplotypeCaller) from input_preprocessed_bams_forHaplotypeCaller.combine(reference_genome_bundle_forHaplotypeCaller)

	output:
	tuple path(gvcf_raw), path(gvcf_raw_index) into raw_gvcfs

	script:
	gvcf_raw = "${bam_preprocessed}".replaceFirst(/\..*bam/, ".g.vcf")
	gvcf_raw_index = "${gvcf_raw}.idx"
	"""
	gatk HaplotypeCaller \
	--java-options "-Xmx14G -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1 -XX:+AggressiveOpts" \
	--verbosity ERROR \
	--max-alternate-alleles 3 \
	--standard-min-confidence-threshold-for-calling 50 \
	-ERC GVCF \
	-R "${reference_genome_fasta_forHaplotypeCaller}" \
	-I "${bam_preprocessed}" \
	-O "${gvcf_raw}"
	"""
}

/*

// GATK GenomicsDBImport ~ import GVCF into data storage system to make data more accessible to tools
process createGenomicsDb_gatk {
	
	input:
	tuple path(gvcf_raw), path(gvcf_raw_index) from raw_gvcfs

	output:
	path genomics_db_workspace into joint_genotyping_db

	script:
	genomics_db_workspace = "genomics_db_workspace"
	"""
	gatk GenomicsDBImport \
	--java-options "-Xmx14G -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1 -XX:+AggressiveOpts -Djava.io.tmpdir=/tmp" \
	--verbosity ERROR \
	-V ${gvcf_raw.collect()} \
	--tmp-dir=/tmp \
	--genomicsdb-workspace-path "${genomics_db_workspace}"
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
	path genomics_db_workspace from joint_genotyping_db
	tuple path(reference_genome_fasta_forJointGenotyping), path(reference_genome_fasta_index_forJointGenotyping), path(reference_genome_fasta_dict_forJointGenotyping) from reference_genome_bundle_forJointGenotyping
	tuple path(gatk_bundle_dbsnp138), path(gatk_bundle_dbsnp138_index) from gatk_reference_bundle_forJointGenotyping

	output:
	tuple path(vcf_joint_genotyped), path(vcf_joint_genotyped_index) into joint_genotyped_vcfs

	script:
	vcf_joint_genotyped = "${params.cohort_name}.vcf"
	vcf_joint_genotyped_index = "${vcf_joint_genotyped}.idx"
	"""
	gatk GenotypeGVCFs \
	--java-options "-Xmx16G -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1 -XX:+AggressiveOpts" \
	--verbosity ERROR \
	--dbsnp "${gatk_bundle_dbsnp138}" \
	--reference "${reference_genome_fasta_forJointGenotyping}" \
	-V gendb://"${genomics_db_workspace}"
	--annotation-group StandardAnnotation \
	--standard-min-confidence-threshold-for-calling 50 \
	-O "${vcf_joint_genotyped}"
	"""
}

// GATK VariantFiltration/SelectVariants ~ hard filter test for excess heterozygosity
// ExcessHet is a phred-scaled p-value. We want a cutoff of anything more extreme
// than a z-score of -4.5 which is a p-value of 3.4e-06, which phred-scaled is 54.69
process excessHeterozygosityHardFilter_gatk {
	publishDir "${params.output_dir}/germline/excessHetHardFilter", mode: 'symlink'

	input:
	tuple path(vcf_joint_genotyped), path(vcf_joint_genotyped_index) from joint_genotyped_vcfs

	output:
	tuple path(vcf_genotyped_filtered), path(vcf_genotyped_filtered_index) into genotyped_filtered_vcfs_forIndelVariantRecalibration, genotyped_filtered_vcfs_forSnvVariantRecalibration

	script:
	intermediate_vcf_genotyped_filtered = "${vcf_joint_genotyped}".replaceFirst(/\.vcf/, ".markedforfilter.vcf")
	vcf_genotyped_filtered = "${vcf_joint_genotyped}".replaceFirst(/\.vcf/, ".filtered.vcf")
	vcf_genotyped_filtered_index = "${vcf_genotyped_filtered}.idx"
	"""
	gatk VariantFiltration \
	--java-options "-Xmx8G -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1 -XX:+AggressiveOpts" \
	--filter-name ExcessHet \
	--filter-expression "ExcessHet > 54.69" \
	-V "${vcf_joint_genotyped}" \
	-O "${intermediate_vcf_genotyped_filtered}"

	gatk SelectVariants \
	--java-options "-Xmx8G -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1 -XX:+AggressiveOpts" \
	--exclude-filtered \
	-V "${intermediate_vcf_genotyped_filtered}" \
	-O "${vcf_genotyped_filtered}"
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

	input:
	tuple path(vcf_genotyped_filtered), path(vcf_genotyped_filtered_index) from genotyped_filtered_vcfs_forIndelVariantRecalibration
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
	--java-options "-Xmx16G -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1 -XX:+AggressiveOpts" \
	--verbosity ERROR \
	--mode INDEL \
	-an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum -an InbreedingCoeff \
	--max-gaussians 4 \
	--trust-all-polymorphic \
	--reference "${reference_genome_fasta_forIndelVariantRecalibration}" \
	--resource:mills,known=false,training=true,truth=true,prior=12 "${gatk_bundle_mills_1000G}" \
	--resource:axiomPoly,known=false,training=true,truth=false,prior=10 "${gatk_bundle_axiom}" \
	--resource:dbsnp,,known=true,training=false,truth=false,prior=2 "${gatk_bundle_dbsnp138_forIndelVariantRecalibration}" \
	-V "${vcf_genotyped_filtered}" \
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

reference_genome_fasta_forSnvVariantRecalibration.combine( reference_genome_fasta_index_forSnvVariantRecalibration )
	.combine( reference_genome_fasta_dict_forSnvVariantRecalibration )
	.set{ reference_genome_bundle_forSnvVariantRecalibration }

// GATK VariantRecalibrator (SNPs) ~ build recalibration model to score SNP variant quality for filtering
process snpVariantRecalibration_gatk {
	publishDir "${params.output_dir}/germline/snpVariantRecal", mode: 'symlink'

	input:


}

*/