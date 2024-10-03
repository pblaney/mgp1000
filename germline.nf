// Myeloma Genome Pipeline 1000
// Comprehensive pipeline for analysis of matched T/N Multiple Myeloma WGS data
// https://github.com/pblaney/mgp1000

// This module of the pipeline is used for germline variant analysis of normal WGS samples.

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
	  nextflow run germline.nf --run_id STR --sample_sheet FILE -profile germline [-bg] [-resume]
      [--input_dir PATH] [--output_dir PATH] [--email STR] [--fastngsadmix STR] [--seq_protocol STR]
      [--cpus INT] [--memory STR] [--queue_size INT] [--executor STR] [--help]

	Mandatory Arguments:
	  --run_id                       STR  Unique identifier for pipeline run
	  --sample_sheet                FILE  CSV file containing the list of samples where the
	                                      first column designates the file name of the normal
	                                      sample, the second column for the file name of the
	                                      matched tumor sample
	  -profile                       STR  Configuration profile to use, must use germline

	Optional Arguments:
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
	                                      and WES for whole-exome, PANEL for MGP targeted panel
	                                      [Default: WGS | Available: WGS, WES, PANEL]
	  --deepvariant                  STR  Indicates whether or not to run DeepVariant workflow
	                                      [Default: on | Available: on, off]
	  --fastngsadmix                 STR  Indicates whether or not to run fastNGSadmix workflow
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

if( params.sample_sheet == null ) exit 1, "The run command issued does not have the '--sample_sheet' parameter set. Please set the '--sample_sheet' parameter to the path of the normal/tumor pair sample sheet CSV."

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

Channel
    .value( file('references/hg38/mgp_panel_calling_regions_v22.hg38.bed') )
    .set{ target_regions_bed }


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
        .into{ input_preprocessed_bams_forDeepVariant;
               input_preprocessed_bams_forDeepVariantPanel;
               input_preprocessed_bams_forAngsd }
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
    params.deepvariant == "on" && params.seq_protocol != "PANEL"

    script:
    sample_id = "${bam_preprocessed}".replaceFirst(/\.final\..*bam/, "")
    deepvariant_germline_vcf_per_chrom = "${sample_id}.deepvariant.germline.snp.indel.${chromosome}.vcf.gz"
    deepvariant_germline_vcf_per_chrom_index = "${deepvariant_germline_vcf_per_chrom}.tbi"
    per_chromosome_bed_file = "${target_bed}".replaceFirst(/\.bed/, ".${chromosome}.bed")
    """
    run_deepvariant \
    --model_type "${params.seq_protocol}" \
    --ref "${ref_genome_fasta}" \
    --reads "${bam_preprocessed}" \
    --regions "${chromosome}" \
    --output_vcf "${sample_id}.deepvariant.raw.germline.snp.indel.${chromosome}.vcf.gz" \
    --output_gvcf "${sample_id}.deepvariant.raw.germline.snp.indel.${chromosome}.gvcf.gz" \
    --num_shards ${task.cpus}

    zcat "${sample_id}.deepvariant.raw.germline.snp.indel.${chromosome}.vcf.gz" \
        | awk '(\$5 !~ ",")' \
        | grep -E '^#|PASS' \
        | bgzip > "${deepvariant_germline_vcf_per_chrom}"

    tabix "${deepvariant_germline_vcf_per_chrom}"
    """
}

// DeepVariant ~ deep learning-based germline variant-calling
process snpAndIndelCallingForPanel_deepvariant {
    tag "${sample_id}"

    input:
    tuple path(bam_preprocessed), path(bam_preprocessed_index) from input_preprocessed_bams_forDeepVariantPanel
    path ref_genome_fasta from ref_genome_fasta_file
    path ref_genome_fasta_index from ref_genome_fasta_index_file
    path ref_genome_fasta_dict from ref_genome_fasta_dict_file
    path target_bed from target_regions_bed

    output:
    tuple val(sample_id), path(deepvariant_germline_panel_vcf), path(deepvariant_germline_panel_vcf_index) into unfiltered_deepvariant_panel_vcf

    when:
    params.deepvariant == "on" && params.seq_protocol == "PANEL"

    script:
    sample_id = "${bam_preprocessed}".replaceFirst(/\.final\..*bam/, "")
    deepvariant_germline_panel_vcf = "${sample_id}.deepvariant.germline.snp.indel.vcf.gz"
    deepvariant_germline_panel_vcf_index = "${deepvariant_germline_panel_vcf}.tbi"
    """
    grep "mutations_and_copynumber" "${target_bed}" | cut -f 1-3 > mgp_panel.bed

    run_deepvariant \
        --model_type WES \
        --ref "${ref_genome_fasta}" \
        --reads "${bam_preprocessed}" \
        --regions mgp_panel.bed \
        --output_vcf "${deepvariant_germline_panel_vcf}" \
        --output_gvcf "${sample_id}.deepvariant.raw.germline.snp.indel.gvcf.gz" \
        --num_shards ${task.cpus}
    """
}

// GATK SortVcfs ~ merge all GVCF files for each sample and sort them
process mergeAndSortVcfs_gatk {
    tag "${sample_id}"
    
    input:
    tuple val(sample_id), path(deepvariant_germline_vcf_per_chrom), path(deepvariant_germline_vcf_per_chrom_index) from unfiltered_germline_vcf.groupTuple()

    output:
    tuple val(sample_id), path(vcf_merged_unfiltered), path(vcf_merged_unfiltered_index) into merged_unfiltered_deepvariant_vcf

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

if( params.seq_protocol == "WGS" | params.seq_protocol == "WES" ) {
    merged_unfiltered_germline_vcf = merged_unfiltered_deepvariant_vcf
} else if( params.seq_protocol == "PANEL" ) {
    merged_unfiltered_germline_vcf = unfiltered_deepvariant_panel_vcf
}

// GATK VariantFiltration ~ hard filter sites based on variant allele frequency and genotype quality
process variantFilter_gatk {
    publishDir "${params.output_dir}/germline", mode: 'copy'
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(vcf_merged_unfiltered), path(vcf_merged_unfiltered_index) from merged_unfiltered_germline_vcf

    output:
    path snp_indel_vcf_merged
    path snp_indel_vcf_merged_index

    when:
    params.deepvariant == "on"

    script:
    snp_indel_vcf_merged = "${sample_id}.deepvariant.germline.snp.indel.vcf.gz"
    snp_indel_vcf_merged_index = "${snp_indel_vcf_merged}.tbi"
    """
    gatk SelectVariants \
    --java-options "-Xmx${task.memory.toGiga()}G -Djava.io.tmpdir=. -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
    --verbosity ERROR \
    --tmp-dir . \
    --exclude-filtered \
    --variant "${vcf_merged_unfiltered}" \
    --output "${snp_indel_vcf_merged}"
    """
}

// ANGSD ~ estimate allele frequencies from genotype likelihoods
process alleleFrequencyEstimation_angsd {
    tag "${sample_id}"

    input:
    tuple path(normal_bam), path(normal_bam_index) from input_preprocessed_bams_forAngsd
    path admix_markers_sites
    path admix_markers_sites_bin
    path admix_markers_sites_index
    //path target_bed from target_regions_bed

    output:
    tuple val(sample_id), path(beagle_genotype_likelihoods) into beagle_input_forFastNgsAdmix

    when:
    params.fastngsadmix == "on"

    script:
    sample_id = "${normal_bam}".replaceFirst(/\.final\..*bam/, "")
    beagle_genotype_likelihoods = "${sample_id}.beagle.gz"

    //if ( params.seq_protocol == "WGS" )
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
    //else ( params.seq_protocol == "WES" | params.seq_protocol == "PANEL" )
    //"""
    //grep -w "${chromosome}" ${target_bed} | grep "mutations_and_copynumber" > "${per_chromosome_bed_file}"
    //"""
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
