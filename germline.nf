// Myeloma Genome Project 1000
// Comprehensive pipeline for analysis of matched T/N Multiple Myeloma WGS data
// https://github.com/pblaney/mgp1000

// This portion of the pipeline is used for germline variant analysis of normal WGS samples.
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

	                          ░█▀▀░█▀▀░█▀▄░█▄█░█░░░▀█▀░█▀█░█▀▀
                              ░█░█░█▀▀░█▀▄░█░█░█░░░░█░░█░█░█▀▀
                              ░▀▀▀░▀▀▀░▀░▀░▀░▀░▀▀▀░▀▀▀░▀░▀░▀▀▀

	Usage:
	  nextflow run germline.nf --run_id STR --sample_sheet FILE -profile germline
	  [-bg] [-resume] [--input_dir PATH] [--output_dir PATH] [--email STR] [--fastngsadmix STR]
	  [--seq_protocol STR] [--cpus INT] [--memory STR] [--queue_size INT] [--executor STR] [--help]

	Mandatory Arguments:
	  --run_id                       STR  Unique identifier for pipeline run
	  --sample_sheet                FILE  CSV file containing the list of samples where the
	                                      first column designates the file name of the normal
	                                      sample, the second column for the file name of the
	                                      matched tumor sample
	  -profile                       STR  Configuration profile to use, must use germline

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
	  --seq_protocol                 STR  Sequencing protocol of the input, WGS for whole-genome
	                                      and WES for whole-exome
	                                      [Default: WGS | Available: WGS, WES]
	  --deepvariant                  STR  Indicates whether or not to run DeepVariant workflow
	                                      [Default: on | Available: on, off]
	  --fastngsadmix                 STR  Indicates whether or not to only run the fastNGSadmix
	                                      workflow
	                                      [Default: on | Available: on, off]
	  --cpus                         INT  Globally set the number of cpus to be allocated
	  --memory                       STR  Globally set the amount of memory to be allocated,
	                                      written as '##.GB' or '##.MB'
	  --queue_size                   INT  Set max number of tasks the pipeline will launch
	                                      [Default: 100]
	  --executor                     STR  Set the job executor for the run
	                                      [Default: slurm | Available: local, slurm, lsf]
	  --help                        FLAG  Prints this message

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
params.deepvariant = "on"
params.fastngsadmix = "on"
//params.fastngsadmix_only = "no"
params.cpus = null
params.memory = null
params.queue_size = 100
params.executor = 'slurm'
params.help = null

// Print help message if requested
if( params.help ) exit 0, helpMessage()

// Print error message if user-defined input/output directories does not exist
if( !file(params.input_dir).exists() ) exit 1, "The user-specified input directory does not exist in filesystem."

// Print error messages if required parameters are not set
if( params.run_id == null ) exit 1, "The run command issued does not have the '--run_id' parameter set. Please set the '--run_id' parameter to a unique identifier for the run."

if( params.sample_sheet == null & params.fastngsadmix_only == "no" ) exit 1, "The run command issued does not have the '--sample_sheet' parameter set. Please set the '--sample_sheet' parameter to the path of the normal/tumor pair sample sheet CSV."

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
	.fromList( ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6',
	            'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
	            'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18',
	            'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY'] )
	.set{ chromosome_list_forDeepVariant }







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


Channel
    .value( file('references/hg38/fastNGSadmix_markers_refPanel.hg38.txt') )
    .set{ admix_markers_ref_population_allele_freqs }

Channel
    .value( file('references/hg38/fastNGSadmix_markers_nInd.hg38.txt') )
    .set{ admix_markers_ref_population_num_of_individuals }

Channel
    .value( file('references/hg38/fastNGSadmix_markers.hg38.sites') )
    .set{ admix_markers_sites }

Channel
    .value( file('references/hg38/fastNGSadmix_markers.hg38.sites.bin') )
    .set{ admix_markers_sites_bin }

Channel
    .value( file('references/hg38/fastNGSadmix_markers.hg38.sites.idx') )
    .set{ admix_markers_sites_index }


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
log.info "        ░█▀▀░█▀▀░█▀▄░█▄█░█░░░▀█▀░█▀█░█▀▀        "
log.info "        ░█░█░█▀▀░█▀▄░█░█░█░░░░█░░█░█░█▀▀        "
log.info "        ░▀▀▀░▀▀▀░▀░▀░▀░▀░▀▀▀░▀▀▀░▀░▀░▀▀▀        "
log.info ''
log.info "~~~ Launch Time ~~~"
log.info ''
log.info ' ${workflowTimestamp}'
log.info ''
log.info "~~~ Input Directory ~~~"
log.info ''
log.info ' ${params.input_dir}'
log.info ''
log.info "~~~ Output Directory ~~~"
log.info ''
log.info ' ${params.output_dir}'
log.info ''
log.info "~~~ Run Report File ~~~"
log.info ''
log.info ' nextflow_report.${params.run_id}.html'
log.info ''
log.info "~~~ Sequencing Protocol ~~~"
log.info ''
log.info " ${params.seq_protocol}"
log.info ''
log.info '################################################'
log.info ''

// Read user provided sample sheet to find Normal sample BAM files
if(params.sample_sheet != null) {
    Channel
        .fromPath( params.sample_sheet )
        .splitCsv( header:true )
        .map{ row -> normal_bam = "${row.normal}"
              		 normal_bam_index = "${row.normal}".replaceFirst(/\.bam$/, "")
              		 return[ file("${params.input_dir}/${normal_bam}"), 
                      		 file("${params.input_dir}/${normal_bam_index}*.bai") ] }
        .unique()
        .into{ input_preprocessed_bams_forAngsd;
               input_preprocessed_bams_forHaplotypeCaller;
               input_preprocessed_bams_forDeepVariant }
}	

// DeepVariant ~ deep learning-based germline variant-calling 
process snpAndIndelCalling_deepvariant {
    tag "${sample_id} C=${chromosome}"

    input:
    tuple path(bam_preprocessed), path(bam_preprocessed_index) from input_preprocessed_bams_forDeepVariant
    path ref_genome_fasta from ref_genome_fasta_file
    path ref_genome_fasta_index from ref_genome_fasta_index_file
    path ref_genome_fasta_dict from ref_genome_fasta_dict_file
    each chromosome from chromosome_list_forDeepVariant

    output:
    tuple val(sample_id), path(deepvariant_germline_vcf_per_chrom), path(deepvariant_germline_vcf_per_chrom_index) into unfiltered_germline_vcf

    when:
    params.deepvariant == "on"

    script:
    sample_id = "${bam_preprocessed}".replaceFirst(/\.final\..*bam/, "")
    deepvariant_germline_vcf_per_chrom = "${sample_id}.deepvariant.germline.snp.indel.${chromosome}.vcf.gz"
    deepvariant_germline_vcf_per_chrom_index = "${deepvariant_germline_vcf_per_chrom}.tbi"
    """
    run_deepvariant \
    --model_type "${params.seq_protocol}" \
    --ref "${ref_genome_fasta}" \
    --reads "${bam_preprocessed}" \
    --regions "${chromosome}" \
    --output_vcf "${sample_id}.deepvariant.raw.germline.snp.indel.${chromosome}.vcf.gz" \
    --num_shards ${task.cpus}

    zcat "${sample_id}.deepvariant.raw.germline.snp.indel.${chromosome}.vcf.gz" \
    | \
    awk '(\$5 !~ ",")' \
    | \
    grep -E '^#|PASS' \
    | \
    gzip > "${deepvariant_germline_vcf_per_chrom}"

    tabix "${deepvariant_germline_vcf_per_chrom}"
    """
}

// GATK SortVcfs ~ merge all GVCF files for each sample and sort them
process mergeAndSortVcfs_gatk {
    tag "${sample_id}"
    
    input:
    tuple val(sample_id), path(deepvariant_germline_vcf_per_chrom), path(deepvariant_germline_vcf_per_chrom_index) from unfiltered_germline_vcf.groupTuple()

    output:
    tuple val(sample_id), path(vcf_merged_unfiltered), path(vcf_merged_unfiltered_index) into merged_unfiltered_germline_vcf

    when:
    params.deepvariant == "on"

    script:
    vcf_merged_unfiltered = "${sample_id}.deepvariant.germline.snp.indel.unfiltered.vcf.gz"
    vcf_merged_unfiltered_index = "${vcf_merged_unfiltered}.tbi"
    """
    gatk SortVcf \
    --java-options "-Xmx${task.memory.toGiga()}G -Djava.io.tmpdir=. -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
    --VERBOSITY ERROR \
    --TMP_DIR . \
    --MAX_RECORDS_IN_RAM 4000000 \
    ${deepvariant_germline_vcf_per_chrom.collect { "--INPUT $it " }.join()} \
    --OUTPUT "${vcf_merged_unfiltered}"
    """
}

// GATK VariantFiltration ~ hard filter sites based on variant allele frequency and genotype quality
process variantFilter_gatk {
    publishDir "${params.output_dir}/germline/${sample_id}", mode: 'copy', pattern: '.{vcf.gz, vcf.gz.tbi}'
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(vcf_merged_unfiltered), path(vcf_merged_unfiltered_index) from merged_unfiltered_germline_vcf

    output:
    path snp_vcf_merged
    path snp_vcf_merged_index
    path indel_vcf_merged
    path indel_vcf_merged_index

    when:
    params.deepvariant == "on"

    script:
    snp_vcf_merged = "${sample_id}.deepvariant.germline.snp.vcf.gz"
    snp_vcf_merged_index = "${snp_vcf_merged}.tbi"
    indel_vcf_merged = "${sample_id}.deepvariant.germline.indel.vcf.gz"
    indel_vcf_merged_index = "${indel_vcf_merged}.tbi"
    """
    gatk VariantFiltration \
    --java-options "-Xmx${task.memory.toGiga()}G -Djava.io.tmpdir=. -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
    --verbosity ERROR \
    --tmp-dir . \
    --filter-name "TopOrBottomQuartileVAF" \
    --filter-expression "(VAF < 0.25 && VAF > 0.0) || VAF > 0.75" \
    --variant "${vcf_merged_unfiltered}" \
    --output "${sample_id}.deepvariant.germline.snp.indel.marked.vcf.gz"

    gatk SelectVariants \
    --java-options "-Xmx${task.memory.toGiga()}G -Djava.io.tmpdir=. -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
    --verbosity ERROR \
    --tmp-dir . \
    --select-type-to-include SNP \
    --exclude-filtered \
    --variant "${sample_id}.deepvariant.germline.snp.indel.marked.vcf.gz" \
    --output "${snp_vcf_merged}"

    gatk SelectVariants \
    --java-options "-Xmx${task.memory.toGiga()}G -Djava.io.tmpdir=. -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
    --verbosity ERROR \
    --tmp-dir . \
    --select-type-to-include INDEL \
    --exclude-filtered \
    --variant "${sample_id}.deepvariant.germline.snp.indel.marked.vcf.gz" \
    --output "${indel_vcf_merged}"
    """
}











/*


// Combine all needed reference FASTA files into one channel for use in SplitIntervals process
reference_genome_fasta_forSplitIntervals.combine( reference_genome_fasta_index_forSplitIntervals )
	.combine( reference_genome_fasta_dict_forSplitIntervals )
	.set{ reference_genome_bundle_forSplitIntervals }

// GATK SplitIntervals ~ divide interval list into files containing equal number of reads for better pipeline performance
process splitIntervalList_gatk {
	
	input:
	tuple path(reference_genome_fasta_forSplitIntervals), path(reference_genome_fasta_index_forSplitIntervals), path(reference_genome_fasta_dict_forSplitIntervals) from reference_genome_bundle_forSplitIntervals
	path gatk_bundle_wgs_interval_list

	output:
	path "splitIntervals/*-split.interval_list" into split_intervals mode flatten

	when:
	params.fastngsadmix_only == "no"

	script:
	"""
	gatk SplitIntervals \
	--java-options "-Xmx${task.memory.toGiga()}G -Djava.io.tmpdir=." \
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
	tuple path(bam_preprocessed), path(bam_preprocessed_index), path(reference_genome_fasta_forHaplotypeCaller), path(reference_genome_fasta_index_forHaplotypeCaller), path(reference_genome_fasta_dict_forHaplotypeCaller), path(interval) from input_bams_and_reference_fasta_forHaplotypeCaller.combine(split_intervals)

	output:
	tuple val(sample_id), path(gvcf_per_interval_raw), path(gvcf_per_interval_raw_index) into raw_gvcfs

	when:
	params.fastngsadmix_only == "no"

	script:
	sample_id = "${bam_preprocessed}".replaceFirst(/\.final\..*bam/, "")
	interval_id = "${interval}".replaceFirst(/-split\.interval_list/, "_int")
	gvcf_per_interval_raw = "${sample_id}.${interval_id}.g.vcf.gz"
	gvcf_per_interval_raw_index = "${gvcf_per_interval_raw}.tbi"
	"""
	gatk HaplotypeCaller \
	--java-options "-Xmx${task.memory.toGiga()}G -Djava.io.tmpdir=." \
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
	tag "${sample_id}"
	
	input:
	tuple val(sample_id), path(gvcf_per_interval_raw), path(gvcf_per_interval_raw_index) from raw_gvcfs.groupTuple()

	output:
	path gvcf_merged_raw into merged_raw_gcvfs
	path gvcf_merged_raw_index into merged_raw_gcvfs_indicies

	when:
	params.fastngsadmix_only == "no"

	script:
	gvcf_merged_raw = "${sample_id}.g.vcf.gz"
	gvcf_merged_raw_index = "${gvcf_merged_raw}.tbi"
	"""
	gatk SortVcf \
	--java-options "-Xmx${task.memory.toGiga()}G -Djava.io.tmpdir=." \
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
	tag "${params.cohort_name}"

	input:
	path gvcf_merged_raw from merged_raw_gcvfs.toList()
	path gvcf_merged_raw_index from merged_raw_gcvfs_indicies.collect()
	tuple path(reference_genome_fasta_forCombineGvcfs), path(reference_genome_fasta_index_forCombineGvcfs), path(reference_genome_fasta_dict_forCombineGvcfs) from reference_genome_bundle_forCombineGvcfs

	output:
	tuple path(gvcf_cohort_combined), path(gvcf_cohort_combined_index) into combined_cohort_gvcf

	when:
	params.fastngsadmix_only == "no"

	script:
	gvcf_cohort_combined = "${params.cohort_name}.g.vcf.gz"
	gvcf_cohort_combined_index = "${gvcf_cohort_combined}.tbi"
	"""
	gatk CombineGVCFs \
	--java-options "-Xmx24576m -Djava.io.tmpdir=." \
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

// GATK GenotypeGVCFs ~ perform joint genotyping
process jointGenotyping_gatk {
	tag "${params.cohort_name}"

	input:
	tuple path(gvcf_cohort_combined), path(gvcf_cohort_combined_index) from combined_cohort_gvcf
	tuple path(reference_genome_fasta_forJointGenotyping), path(reference_genome_fasta_index_forJointGenotyping), path(reference_genome_fasta_dict_forJointGenotyping) from reference_genome_bundle_forJointGenotyping
	tuple path(gatk_bundle_dbsnp138), path(gatk_bundle_dbsnp138_index) from gatk_reference_bundle_forJointGenotyping

	output:
	tuple path(vcf_joint_genotyped), path(vcf_joint_genotyped_index) into joint_genotyped_vcfs

	when:
	params.fastngsadmix_only == "no"

	script:
	vcf_joint_genotyped = "${params.cohort_name}.vcf.gz"
	vcf_joint_genotyped_index = "${vcf_joint_genotyped}.tbi"
	"""
	gatk GenotypeGVCFs \
	--java-options "-Xmx24576m -Djava.io.tmpdir=." \
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

// GATK VariantFiltration ~ hard filter test for excess heterozygosity
// ExcessHet is a phred-scaled p-value. We want a cutoff of anything more extreme
// than a z-score of -4.5 which is a p-value of 3.4e-06, which phred-scaled is 54.69
process excessHeterozygosityHardFilter_gatk {
	tag "${params.cohort_name}"

	input:
	tuple path(vcf_joint_genotyped), path(vcf_joint_genotyped_index) from joint_genotyped_vcfs

	output:
	tuple path(vcf_hard_filtered), path(vcf_hard_filtered_index) into hard_filtered_vcfs_forIndelVariantRecalibration, hard_filtered_vcfs_forSnpVariantRecalibration, hard_filtered_vcfs_forApplyVqsr

	when:
	params.fastngsadmix_only == "no"

	script:
	vcf_hard_filtered_marked = "${vcf_joint_genotyped}".replaceFirst(/\.vcf\.gz/, ".filtermarked.vcf.gz")
	vcf_hard_filtered = "${vcf_joint_genotyped}".replaceFirst(/\.vcf\.gz/, ".filtered.vcf.gz")
	vcf_hard_filtered_index = "${vcf_hard_filtered}.tbi"
	"""
	gatk VariantFiltration \
	--java-options "-Xmx${task.memory.toGiga()}G -Djava.io.tmpdir=." \
	--verbosity ERROR \
	--tmp-dir . \
	--filter-name ExcessHet \
	--filter-expression "ExcessHet > 54.69" \
	--variant "${vcf_joint_genotyped}" \
	--output "${vcf_hard_filtered_marked}"

	gatk SelectVariants \
	--java-options "-Xmx${task.memory.toGiga()}G -Djava.io.tmpdir=." \
	--verbosity ERROR \
	--tmp-dir . \
	--exclude-filtered \
	--variant "${vcf_hard_filtered_marked}" \
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
	tag "${params.cohort_name}"

	input:
	tuple path(vcf_hard_filtered), path(vcf_hard_filtered_index) from hard_filtered_vcfs_forIndelVariantRecalibration
	tuple path(gatk_bundle_mills_1000G), path(gatk_bundle_mills_1000G_index), path(gatk_bundle_axiom), path(gatk_bundle_axiom_index), path(gatk_bundle_dbsnp138_forIndelVariantRecalibration), path(gatk_bundle_dbsnp138_index_forIndelVariantRecalibration) from gatk_reference_bundle_forIndelVariantRecalibration
	tuple path(reference_genome_fasta_forIndelVariantRecalibration), path(reference_genome_fasta_index_forIndelVariantRecalibration), path(reference_genome_fasta_dict_forIndelVariantRecalibration) from reference_genome_bundle_forIndelVariantRecalibration

	output:
	tuple path(indel_vqsr_table), path(indel_vqsr_table_index), path(indel_vqsr_tranches) into indel_vqsr_files

	when:
	params.fastngsadmix_only == "no"

	script:
	indel_vqsr_table = "${params.cohort_name}.indel.recaldata.table"
	indel_vqsr_table_index = "${indel_vqsr_table}.idx"
	indel_vqsr_tranches = "${params.cohort_name}.indel.tranches"
	"""
	gatk VariantRecalibrator \
	--java-options "-Xmx${task.memory.toGiga()}G -Djava.io.tmpdir=." \
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

// Combine all needed GATK bundle files and reference FASTA files into one channel for use in GATK SNP VariantRecalibratior process
gatk_bundle_hapmap.combine( gatk_bundle_hapmap_index )
	.combine( gatk_bundle_1000G_omni )
	.combine( gatk_bundle_1000G_omni_index )
	.combine( gatk_bundle_1000G_snps )
	.combine( gatk_bundle_1000G_snps_index)
	.combine( gatk_bundle_dbsnp138_forSnpVariantRecalibration )
	.combine( gatk_bundle_dbsnp138_index_forSnpVariantRecalibration )
	.set{ gatk_reference_bundle_forSnpVariantRecalibration }

reference_genome_fasta_forSnpVariantRecalibration.combine( reference_genome_fasta_index_forSnpVariantRecalibration )
	.combine( reference_genome_fasta_dict_forSnpVariantRecalibration )
	.set{ reference_genome_bundle_forSnpVariantRecalibration }

// GATK VariantRecalibrator (SNPs) ~ build recalibration model to score SNP variant quality for filtering
process snpVariantRecalibration_gatk {
	tag "${params.cohort_name}"

	input:
	tuple path(vcf_hard_filtered), path(vcf_hard_filtered_index) from hard_filtered_vcfs_forSnpVariantRecalibration
	tuple path(gatk_bundle_hapmap), path(gatk_bundle_hapmap_index), path(gatk_bundle_1000G_omni), path(gatk_bundle_1000G_omni_index), path(gatk_bundle_1000G_snps), path(gatk_bundle_1000G_snps_index), path(gatk_bundle_dbsnp138_forSnpVariantRecalibration), path(gatk_bundle_dbsnp138_index_forSnpVariantRecalibration) from gatk_reference_bundle_forSnpVariantRecalibration
	tuple path(reference_genome_fasta_forSnpVariantRecalibration), path(reference_genome_fasta_index_forSnpVariantRecalibration), path(reference_genome_fasta_dict_forSnpVariantRecalibration) from reference_genome_bundle_forSnpVariantRecalibration

	output:
	tuple path(snp_vqsr_table), path(snp_vqsr_table_index), path(snp_vqsr_tranches) into snp_vqsr_files

	when:
	params.fastngsadmix_only == "no"

	script:
	snp_vqsr_table = "${params.cohort_name}.snp.recaldata.table"
	snp_vqsr_table_index = "${snp_vqsr_table}.idx"
	snp_vqsr_tranches = "${params.cohort_name}.snp.tranches"	
	"""
	gatk VariantRecalibrator \
	--java-options "-Xmx${task.memory.toGiga()}G -Djava.io.tmpdir=." \
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
	tag "${params.cohort_name}"

	input:
	tuple path(vcf_hard_filtered), path(vcf_hard_filtered_index) from hard_filtered_vcfs_forApplyVqsr
	tuple path(indel_vqsr_table), path(indel_vqsr_table_index), path(indel_vqsr_tranches) from indel_vqsr_files
	tuple path(snp_vqsr_table), path(snp_vqsr_table_index), path(snp_vqsr_tranches) from snp_vqsr_files

	output:
	tuple path(final_vqsr_germline_vcf), path(final_vqsr_germline_vcf_index) into vqsr_germline_vcfs

	when:
	params.fastngsadmix_only == "no"

	script:
	intermediate_vqsr_germline_vcf = "${vcf_hard_filtered}".replaceFirst(/\.filtered\.vcf\.gz/, ".intermediate.vqsr.vcf.gz")
	final_vqsr_germline_vcf = "${vcf_hard_filtered}".replaceFirst(/\.filtered\.vcf\.gz/, ".final.vqsr.vcf.gz")
	final_vqsr_germline_vcf_index = "${final_vqsr_germline_vcf}.tbi"
	"""
	gatk ApplyVQSR \
	--java-options "-Xmx${task.memory.toGiga()}G -Djava.io.tmpdir=." \
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
	--java-options "-Xmx${task.memory.toGiga()}G -Djava.io.tmpdir=." \
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

// Combine all needed reference FASTA files into one channel for use in BCFtools Norm process
reference_genome_fasta_forSplitAndNorm.combine( reference_genome_fasta_index_forSplitAndNorm )
	.combine( reference_genome_fasta_dict_forSplitAndNorm )
	.set{ reference_genome_bundle_forSplitAndNorm }

// BCFtools Norm ~ split multiallelic sites into multiple rows then left-align and normalize indels
process splitMultiallelicAndLeftNormalizeVcf_bcftools {
	tag "${params.cohort_name}"

	input:
	tuple path(final_vqsr_germline_vcf), path(final_vqsr_germline_vcf_index) from vqsr_germline_vcfs
	tuple path(reference_genome_fasta_forSplitAndNorm), path(reference_genome_fasta_index_forSplitAndNorm), path(reference_genome_fasta_dict_forSplitAndNorm) from reference_genome_bundle_forSplitAndNorm

	output:
	tuple path(final_germline_vcf), path(final_germline_vcf_index) into final_germline_vcf_forAnnotation, final_germline_vcf_forAdmixture
	path multiallelics_stats
	path realign_normalize_stats

	when:
	params.fastngsadmix_only == "no"

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
	--output "${final_germline_vcf}"

	tabix "${final_germline_vcf}"
	"""
}




*/



// ANGSD ~ estimate allele frequencies from genotype likelihoods
process alleleFrequencyEstimation_angsd {
    tag "${sample_id}"

    input:
    tuple path(normal_bam), path(normal_bam_index) from input_preprocessed_bams_forAngsd
    path admix_markers_sites
    path admix_markers_sites_bin
    path admix_markers_sites_index

    output:
    tuple val(sample_id), path(beagle_genotype_likelihoods) into beagle_input_forFastNgsAdmix

    when:
    params.fastngsadmix == "on"

    script:
    sample_id = "${normal_bam}".replaceFirst(/\.final\..*bam/, "")
    beagle_genotype_likelihoods = "${sample_id}.beagle.gz"
    """
    angsd \
    -nThreads "${task.cpus}" \
    -i "${normal_bam}" \
    -sites "${admix_markers_sites}" \
    -GL 2 \
    -doGlf 2 \
    -doMajorMinor 3 \
    -minMapQ 30 \
    -minQ 20 \
    -doDepth 1 \
    -doCounts 1 \
    -out "${sample_id}"
    """
}

// fastNGSadmix ~ estimation of sample ancestry using autosomal SNP genotype data in a supervised fashion
process ancestryEstimation_fastngsadmix {
    publishDir "${params.output_dir}/germline/${sample_id}", mode: 'copy'
    tag "${sample_id}"
    
    input:
    tuple val(sample_id), path(beagle_genotype_likelihoods) from beagle_input_forFastNgsAdmix
    path admix_markers_ref_population_allele_freqs
    path admix_markers_ref_population_num_of_individuals

    output:
    path sample_admixture_estimations
    path fastngsadmix_log

    when:
    params.fastngsadmix == "on"

    script:
    sample_admixture_estimations = "${sample_id}.fastngsadmix.23.qopt"
    fastngsadmix_log = "${sample_id}.fastngsadmix.23.log"
    """
    fastNGSadmix \
    -likes "${beagle_genotype_likelihoods}" \
    -fname "${admix_markers_ref_population_allele_freqs}" \
    -Nname "${admix_markers_ref_population_num_of_individuals}" \
    -out "${sample_id}.fastngsadmix.23" \
    -boot 100 \
    -whichPops all
    """
}
