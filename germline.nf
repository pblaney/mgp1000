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

// Set channels for reference files
Channel
	.fromPath( 'references/hg38/bwa/genome.fa' )
	.into{ reference_genome_fasta_forHaplotypeCaller;
	       reference_genome_fasta_forCNNScoreVariants }

Channel
	.fromPath( 'references/hg38/bwa/genome.fa.fai' )
	.into{ reference_genome_fasta_index_forHaplotypeCaller;
		   reference_genome_fasta_index_forCNNScoreVariants }

Channel
	.fromPath( 'references/hg38/bwa/genome.dict' )
	.into{ reference_genome_fasta_dict_forHaplotypeCaller;
	       reference_genome_fasta_dict_forCNNScoreVariants }

Channel
	.fromPath( 'references/hg38/gatkBundle/1000G_Omni2.5.hg38.vcf.gz' )
	.set{ gatk_bundle_1000G_omni }

Channel
	.fromPath( 'references/hg38/gatkBundle/1000G_Omni2.5.hg38.vcf.gz.tbi' )
	.set{ gatk_bundle_1000G_omni_index }

Channel
	.fromPath( 'references/hg38/gatkBundle/Hapmap_3.3.hg38.vcf.gz' )
	.set{ gatk_bundle_hapmap }

Channel
	.fromPath( 'references/hg38/gatkBundle/Hapmap_3.3.hg38.vcf.gz.tbi' )
	.set{ gatk_bundle_hapmap_index }


// #################################################### \\
// ~~~~~~~~~~~~~~~~ PIPELINE PROCESSES ~~~~~~~~~~~~~~~~ \\

// Set all channels for input BAM files
Channel
	.fromPath( 'input/preprocessedBams/*.bam' )
	.into{ input_preprocessed_bams_forTelomereLengthEstimation;
	       input_preprocessed_bams_forHaplotypeCaller;
	       input_preprocessed_bams_forCNNScoreVariants }

// TelSeq ~ estimate telomere length of sample
process telomereLengthEstimation_telseq {
	publishDir "${params.output_dir}/germline/telomereLength", mode: 'move'

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

	input:
	path bam_preprocessed from input_preprocessed_bams_forHaplotypeCaller
	tuple path(reference_genome_fasta_forHaplotypeCaller), path(reference_genome_fasta_index_forHaplotypeCaller), path(reference_genome_fasta_dict_forHaplotypeCaller) from reference_genome_bundle_forHaplotypeCaller

	output:
	path gvcf_raw into raw_gvcfs

	script:
	gvcf_raw = "${bam_preprocessed}".replaceFirst(/\..*bam/, ".raw.gvcf.gz")
	"""
	gatk HaplotypeCaller \
	--java-options "-Xmx16G -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1 -XX:+AggressiveOpts" \
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


// Combine all needed reference FASTA files into one channel for use in GATK CNNScoreVariants process
reference_genome_fasta_forCNNScoreVariants.combine( reference_genome_fasta_index_forCNNScoreVariants )
	.combine( reference_genome_fasta_dict_forCNNScoreVariants )
	.set{ reference_genome_bundle_forCNNScoreVariants }

// GATK CNNScoreVariants ~ apply a convolutional neural network to score annotated variants for filtering
process cnnScoreVariants_gatk {
	publishDir "${params.output_dir}/germline/cnnScoreVariants", mode: 'symlink'

	input:
	path bam_preprocessed from input_preprocessed_bams_forCNNScoreVariants
	tuple path(reference_genome_fasta_forCNNScoreVariants), path(reference_genome_fasta_index_forCNNScoreVariants), path(reference_genome_fasta_dict_forCNNScoreVariants) from reference_genome_bundle_forCNNScoreVariants
	path vcf_raw from raw_vcfs

	output:
	path vcf_annotated into annotated_vcfs

	script:
	vcf_annotated = "${vcf_raw}".replaceFirst(/\.raw\.vcf\.gz/, ".annotated.vcf.gz")
	"""
	gatk CNNScoreVariants \
	--java-options "-Xmx16G -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1 -XX:+AggressiveOpts" \
	--verbosity ERROR \
	--tensor-type read_tensor \
	-I "${bam_preprocessed}" \
	-R "${reference_genome_fasta_forCNNScoreVariants}" \
	-V "${vcf_raw}" \
	-O "${vcf_annotated}"
	"""
}

// Combine GATK bundle files into one channel for use in GATK FilterVariantTranches process
gatk_bundle_1000G_omni.combine( gatk_bundle_1000G_omni_index )
	.combine( gatk_bundle_hapmap )
	.combine( gatk_bundle_hapmap_index)
	.set{ gatk_reference_bundle }

// GATK FilterVariantTranches ~ apply tranche filtering to VCF based on CNN scores
process trancheFilterVariants_gatk {
	publishDir "${params.output_dir}/germline/trancheFilterVariants", mode: 'symlink'

	input:
	path vcf_annotated from annotated_vcfs
	tuple path(gatk_bundle_1000G_omni), path(gatk_bundle_1000G_omni_index), path(gatk_bundle_hapmap), path(gatk_bundle_hapmap_index) from gatk_reference_bundle

	output:
	path vcf_filtered

	script:
	vcf_filtered = "${vcf_annotated}".replaceFirst(/\.annotated\.vcf\.gz/, ".filtered.vcf.gz")
	"""
	gatk FilterVariantTranches \
	--java-options "-Xmx16G -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1 -XX:+AggressiveOpts" \
	--verbosity ERROR \
	--info-key CNN_2D \
	-V "${vcf_annotated}" \
	--resource "${gatk_bundle_1000G_omni}" \
	--resource "${gatk_bundle_hapmap}" \
	--snp-tranche 99.9 \
	--indel-tranche 99.6 \
	--invalidate-previous-filters \
	-O "${vcf_filtered}"
	"""
}

*/