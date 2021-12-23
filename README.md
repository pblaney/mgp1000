# Myeloma Genome Project 1000
Comprehensive bioinformatics pipeline for the large-scale collaborative analysis of Multiple Myeloma genomes in an effort to delineate the broad spectrum of somatic events

## Pipeline Overview
In order to analyze over one thousand matched tumor/normal whole-genome samples across multiple data centers in a consistent manner, a pipeline was created that leverages the workflow management, portability, and reproducibility of [Nextflow](http://www.nextflow.io/) in conjuction with [Singularity](https://sylabs.io/docs/).

The entire pipeline is divided into 3 steps: Preprocessing, Germline Variant Analysis, and Somatic Variant Analysis
This compartmentalizes the workflow and provides significant completion checkpoints which is effective for large-scale batch processing. 

<img src="https://github.com/pblaney/mgp1000/blob/master/MGP1000Pipeline.png" width="800">

## Deploying the Pipeline
The pipeline was developed to be run on various HPCs without concern of environment incompatibilities, version issues, or missing dependencies. None of the commands require admin access or `sudo` to be completed. However, there are a few assumptions regarding initial setup of the pipeline but the required software should be readily available on nearly all HPC systems.
* Git
* GNU Utilities
* Java 8 (or later)
* Singularity (validated on v3.1, v3.5.2, v3.7.1 other versions will be tested)


## Clone GitHub Repository
The first step in the deployment process is to clone the MGP1000 GitHub repository to a location on your HPC that is large enough to hold the input/output data and has access to the job scheduling software, such as SLURM or SGE.
```
git clone https://github.com/pblaney/mgp1000.git
```

### Installing Git LFS
In an effort to containerize the pipeline further, all the necessary reference files and Singularity container images are stored in the GitHub repository using their complementary [Large File Storage (LFS)](https://git-lfs.github.com) extension. This requires a simple installation of the binary executible file at a location on your `$PATH`. The extension pairs seemlessly with Git to download all files while cloning the repository.

**NOTE: Many HPC environments may already have this dependency installed, if so this section can be skipped.**

If required, a `make` command will complete the installation of Linux AMD64 binary executible git-lfs file (v3.0.2). Other binary files available [here](https://github.com/git-lfs/git-lfs/releases)


```
make install-gitlfs-linuxamd64
```
Move the `git-lfs` binary to a location on `$PATH`
```
mv git-lfs $HOME/bin
```
Use `git-lfs` to retrieve the rest of the GitHub repository's files
```
git-lfs pull
```

### Reference Data
To facilitate ease of use, reproducibility, and consistency between all users of the pipeline, all required reference data has been provided within the `references/hg38/` directory. Detailed provenance of each file per tool is included in the pipline [Wiki](https://github.com/pblaney/mgp1000/wiki) for full traceability.

### Containers
For the same reasons as with the reference data, the Singularity image files needed for each tool's container is provided within the `containers/` directory. All containers were originally developed with Docker and all tags can be found on the associated [DockerHub](https://hub.docker.com/r/patrickblaneynyu/mgp1000)

## Install Nextflow
This series of commands will first check if Java is available to the base pipeline environment and then install Nextflow.

**NOTE: Many HPC environments may already have this dependency installed, if so this section can be skipped.**
```
java -version
```
```
make install-nextflow
```
For ease of use, ove the binary executible `nextflow` file to same directory as git-lfs
```
mv nextflow $HOME/bin
```


## Prepare the Pipeline for Usage
Due to size, certain reference files are GNU zipped so the `make prep-pipeline` command must be run to prepare them for use in the pipeline. Additionally, an `input` directory is created for staging all input BAM or FASTQ files, `preprocessedBams` subdirectory for BAMs that have undergone preprocessing and are ready for Germline/Somatic Variant Analysis steps, and a `logs` directory to store Nextflow output log files for each run.
```
make prep-pipeline
```


## Stage Input BAM or FASTQ Files
By default, all input files are handled out of the `input` and `input/preprocessedBams` directories for the Preprocessing and Germline/Somatic Variant Analysis steps, respectively. However, each step in the pipeline includes an option (`--input_dir`) for the user to define the input directory. Additionally, the pipeline will follow symbolic links for input files so there is no need to move files for staging. Given the possible size of the input data, the samples may have to be processed in batches. Additionally, the pipeline is designed to process batches of identical format, i.e. all BAMs or all FASTQs.

**NOTE: A key assumption is that any input FASTQs use an 'R1/R2' naming convention to designate paired-end read files. Check the `testSample` directory to see examples of FASTQ naming conventions that are accepted. It is recommended that these be used as a sanity check of the pipeline if deploying for the first time.**

Example of staging input data files with symbolic link
```
ln -s /absolute/path/to/unprocessed/samples/directory/*.fastq.gz input/
```


## Run the Preprocessing Step of the Pipeline
The Preprocessing step of the pipeline will be started with one command that will handle linking each individual process in the pipeline to the next. A key advantage of using Nextflow within an HPC environment is that will also perform all the job scheduling/submitting given the correct configuration with the user's [executor](https://www.nextflow.io/docs/latest/executor.html).

There are two methods for running each step in the pipeline: [directly from the current environment](#Direct Submission Example) or submitted to job scheduler as a whole. Typically, it is regarded as best practice in an HPC setting to submit the job as a whole but both are equally handled and output is the same. Examples of both are given below.

**NOTE: The pipeline is currently configured to run with SLURM as the executor. If the user's HPC uses an alternative scheduler please reach out for assistance with adjustments to the configuration to accommodate this, contact information at end of README.**

### Direct Submission Example
First, the user will need to load the necessary software to run the pipeline step to the environment. At most, this will require Java, Nextflow, and Singularity. This is a user-specific step so the commands may be different depending on the user's HPC configuration.
```
$ module load java/1.8 nextflow/21.04.3 singularity/3.7.1
```
Basic example of direct submission after loading in required modules for execution.
```
$ nextflow run preprocessing.nf -bg --run_id batch1 --input_format fastq --email someperson@gmail.com -profile preprocessing
```
Here is the full help message for the Preprocessing step.
```
$ nextflow run preprocessing.nf --help
...
...
...
################################################

Usage Example:

	nextflow run preprocessing.nf -bg -resume --run_id batch1 --input_format fastq --email someperson@gmail.com -profile preprocessing 

Mandatory Arguments:
	--run_id                       [str]  Unique identifier for pipeline run
	--input_format                 [str]  Format of input files
	                                      Available: fastq, bam
	-profile                       [str]  Configuration profile to use, each profile described in nextflow.config file
	                                      Available: preprocessing, germline, somatic

Main Options:
	-bg                           [flag]  Runs the pipeline processes in the background, this option should be included if deploying
	                                      pipeline with real data set so processes will not be cut if user disconnects from deployment
	                                      environment
	-resume                       [flag]  Successfully completed tasks are cached so that if the pipeline stops prematurely the
	                                      previously completed tasks are skipped while maintaining their output
	--input_dir                    [str]  Directory that holds BAMs and associated index files
	                                      Default: input/
	--output_dir                   [str]  Directory that will hold all output files from the somatic variant analysis
	                                      Default: output/
	--email                        [str]  Email address to send workflow completion/stoppage notification
	--skip_to_qc                   [str]  Skips directly to final Preprocessing QC step, can only be used in conjunction with bam as the input_format,
	                                      should only be used for extreme coverage BAMs that have been previously aligned with BWA MEM to the hg38
	                                      reference genome and have adequate provenance to reflect this
	                                      Available: yes, no
	                                      Default: no
	--cpus                         [int]  Globally set the number of cpus to be allocated for all processes
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

################################################
```

### Preprocessing Output Description
This step of the pipeline generates various per tool QC metrics that are useful in determining samples for best use in downstream analyses. By default, all output files are stored into process-specific subdirectories within the `output/` directory. However, each step in the pipeline includes an option (`--output_dir`) for the user to define the base output directory.

Here is a snapshot of the expected subdirectories within the `output/preprocessing` base directory after a successful run of the Preprocessing step:

| Subdirectory | Output Files | Description of Files |
| --- | --- | --- |
| `trimLogs` | `*.trim.log` | number of reads before and after trimming for quality |
| `fastqc` | `*_fastqc.[html / zip]` | in-depth quality evaluation on a per base and per sequence manner | 
| `alignmentFlagstats` | `*.alignment.flagstat.log` | initial number of reads of various alignment designation |
| `markdupFlagstats` | `*.markdup.[log/flagstat.log]` | number of detected duplicate reads, number of reads after deduplication |
| `finalPreprocessedBams` | `*.final.[bam/bai]` | final preprocessed BAM and index for downstream analysis |
| `coverageMetrics` | `*.coverage.metrics.txt` | genome-wide coverage metrics of final BAM |
| `gcBiasMetrics` | `*.gcbias.[metrics.txt/metrics.pdf/summary.txt]` | genome-wide GC bias metrics of final BAM |

Upon completion of the Preprocessing run, there is a `make preprocessing-completeion` command that is useful for collecting the run-related logs.
```
make preprocessing-completion
```


## Run the Germline Variant Analysis Step of the Pipeline
The most important component of this step of the pipeline is the user-provided sample sheet CSV. This file includes two comma-separated columns: filename of normal sample BAMs and filename of corresponding paired tumor sample BAMs. An example of this is provided in `samplesheet.csv` within the `testSamples` directory. The sample sheet file should typically be within the main `mgp1000` directory.

### Note on Parameters
There are two parameters that will prepare necessary reference files as part of this step of the pipeline, `--vep_ref_cached` and `--ref_vcf_concatenated`. These parameters must be be set to `no` for the first run of the Germline Variant Analysis step of the pipeline.  
```
nextflow run germline.nf --help
...
...
...
################################################

Usage Example:

	nextflow run germline.nf -bg -resume --run_id batch1 --sample_sheet samplesheet.csv --cohort_name wgs_set --email someperson@gmail.com --vep_ref_cached no --ref_vcf_concatenated no -profile germline 

Mandatory Arguments:
   	--run_id                       [str]  Unique identifier for pipeline run
   	--sample_sheet                 [str]  CSV file containing the list of samples where the first column designates the file name of the
   	                                      normal sample, the second column for the file name of the matched tumor sample, example of the
   	                                      format for this file is in the testSamples directory
   	--cohort_name                  [str]  A user defined collective name of the group of samples being run through this step of the
   	                                      pipeline, this will be used as the name of the final output multi-sample GVCF
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
	--vep_ref_cached               [str]  Indicates whether or not the VEP reference files used for annotation have been downloaded/cached
	                                      locally, this will be done in a process of the pipeline if it has not, this does not need to be
	                                      done for every separate run after the first
	                                      Available: yes, no
	                                      Default: yes
	--ref_vcf_concatenated         [str]  Indicates whether or not the 1000 Genomes Project reference VCF used for ADMIXTURE analysis has
	                                      been concatenated, this will be done in a process of the pipeline if it has not, this does not
	                                      need to be done for every separate run after the first
	                                      Available: yes, no
	                                      Default: yes
	--cpus                         [int]  Globally set the number of cpus to be allocated for all processes
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

################################################
```

### Germline Variant Analysis Output Description
This step of the pipeline generates a per-cohort joint genotyped VCF and ADMIXTURE estimation of individual ancestries in the context of the 26 populations outlined in the 1000 Genomes Project. By default, all output files are stored into a subdirectory which is named based on the user-defined `--cohort_name` parameter within the `output/` directory. However, each step in the pipeline includes an option (`--output_dir`) for the user to define the base output directory.

Here is a snapshot of the expected subdirectories within the `output/germline/[cohort_name]` base directory after a successful run of the Germline Variant Analysis step:

| Germline Output File | Description of File |
| --- | --- |
| `*.germline.annotated.vcf.gz` | filtered and annotated germline SNP VCF |
| `*.germline.vep.summary.html` | annotation summary HTML file for SNPs |
| `*.hardfiltered.refmerged.stats.txt` | number of SNP sites filtered out before use in ADMIXTURE analysis |
| `*.maf.gt.filtered.refmerged.stats.txt` | number of SNP sites filtered out based on MAF > 0.05 and missing genotypes |
| `*.pruned.maf.gt.filtered.refmerged.stats.txt` | number of SNP sites filtered out based on linkage disequilibrium |
| `*.pruned.maf.gt.filtered.refmerged.pop` | population file used for supervised analysis |
| `*.pruned.maf.gt.filtered.refmerged.26.Q` | ADMIXTURE ancestry fractions |
| `*.pruned.maf.gt.filtered.refmerged.26.P` | ADMIXTURE population allele frequencies |
| `*.pruned.maf.gt.filtered.refmerged.26.Q_se` | standard error of ADMIXTURE ancestry fractions |

Upon completion of the Germline Variant Analysis, there is a `make germline-completeion` command that is useful for collecting the run-related logs.
```
make germline-completion
```


## Run the Somatic Variant Analysis Step of the Pipeline
This step uses the same user-provided sample sheet CSV as the Germline Variant Analysis step.

### Note on Parameters
There are two parameters that will prepare necessary reference files as part of this step of the pipeline, `--vep_ref_cached` and `--mutect_ref_vcf_concatenated`. These parameters must be be set to `no` for the first run of the Somatic Variant Analysis step of the pipeline. By default, all tools in this step will be used with the standard command in the usage example. The consensus output per variant type expects all tools to be included in the run for consensus processes to be run.
```
nextflow run somatic.nf --help
...
...
...
################################################

Usage Example:

	nextflow run somatic.nf -bg -resume --run_id batch1 --sample_sheet samplesheet.csv --email someperson@gmail.com --mutect_ref_vcf_concatenated no --vep_ref_cached no -profile somatic 

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
	--vep_ref_cached               [str]  Indicates whether or not the VEP reference files used for annotation have been downloaded/cached
	                                      locally, this will be done in a process of the pipeline if it has not, this does not need to be
	                                      done for every separate run after the first
	                                      Available: yes, no
	                                      Default: yes
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
	--ascatngs                     [str]  Indicates whether or not to use this tool
	                                      Available: off, on
	                                      Default: on
	--ascatngs_ploidy              [int]  Manually set the ploidy value to be used for the ascatNgs algorithm, this is best used when
	                                      re-running the analysis for specific samples with significant outlier output identified via the
	                                      sunrise plot, setting this also requires setting --ascatngs_purity
	                                      Available: 2, 3, 4, etc.
	                                      Default: estimated internally by tool
	--ascatngs_purity            [float]  Manually set the purity value to be used for the ascatNgs algorithm, this is best used when
	                                      re-running the analysis for specific samples with significant outlier output identified via the
	                                      sunrise plot, setting this also requires setting --ascatngs_ploidy
	                                      Available: 0.2, 0.5, 0.8, etc.
	                                      Default: estimated internally by tool
	--controlfreec                 [str]  Indicates whether or not to use this tool
	                                      Available: off, on
	                                      Default: on
	--sclust                       [str]  Indicates whether or not to use this tool
	                                      Available: off, on
	                                      Default: on
	--sclust_lambda                [str]  Manually set the degree of smoothing for clustering mutations, increasing the value should resolve
	                                      issues with QP iterations related errors
	                                      Available: 1e-6, 1e-5
	                                      Default: 1e-7
	--manta                        [str]  Indicates whether or not to use this tool
	                                      Available: off, on
	                                      Default: on
	--svaba                        [str]  Indicates whether or not to use this tool
	                                      Available: off, on
	                                      Default: on
	--delly                        [str]  Indicates whether or not to use this tool
	                                      Available: off, on
	                                      Default: on

################################################
```

### Somatic Variant Analysis Output Description
This step of the pipeline generates per tumor-normal pair consensus calls for SNVs, InDels, CNVs, and SVs, capture telomere length and composition, and aggregate metadata information on tumor-normal concordance, contamination, purity, ploidy, and subclonal populations. Each tool used has its native output kept within a self-named subdirectory while the final consensus output files per tumor-normal pair are funneled into the `consensus` subdirectory. By default, all output files are stored into process-specific subdirectories within the `output/` directory. However, each step in the pipeline includes an option (`--output_dir`) for the user to define the base output directory.

Here is a snapshot of the final `output/somatic/consensus/[tumor_normal_id]` directory after a successful run of the Somatic Variant Analysis step:

| Somatic Consensus Output File | Description of File |
| --- | --- |
| `*.hq.consensus.somatic.snv.annotated.vcf.gz` | filtered and annotated consensus SNV VCF |
| `*.hq.consensus.somatic.snv.vep.summary.html` | annotation summary HTML file for SNVs |
| `*.hq.consensus.somatic.indel.annotated.vcf.gz` | filtered and annotated consensus InDel VCF |
| `*.hq.consensus.somatic.indel.vep.summary.html` | annotation summary HTML file for InDels |
| `*.consensus.somatic.cnv.alleles.merged.bed` | per segment consensus CNV BED |
| `*.consensus.somatic.cnv.subclonal.txt` | aggregated subclonal population estimates from Control-FREEC and Sclust |
| `*.consensus.somatic.sv.vcf` | consensus SV VCF |
| `*.consensus.somatic.sv.bedpe` | consensus SV calls in simplified BEDPE format |
| `*.consensus.somatic.metadata.txt` | aggregated metadata from alleleCount, Conpair, Mutect2, ascatNGS, Control-FREEC, and Sclust |

Additional per-tool subdirectories included in the base `output/somatic` output directory:

| Subdirectory | Description of Files |
| --- | --- |
| `ascatNgs` | native output of ascatNGS somatic analysis workflow |
| `conpair` | native output of Conpair somatic analysis workflow |
| `controlFreec` | native output of Control-FREEC somatic analysis workflow |
| `delly` | native output of DELLY2 somatic analysis workflow |
| `manta` | native output of Manta somatic analysis workflow |
| `mutect` | native output of Mutect2 somatic analysis workflow |
| `sclust` | native output of Sclust somatic analysis workflow |
| `sexOfSamples` | sample sex estimation using alleleCount |
| `strelka` | native output of Strelka2 somatic analysis workflow |
| `svaba` | native output of SvABA somatic analysis workflow |
| `telomereHunter` | native output of TelomereHunter somatic analysis workflow |
| `varscan` | native output of VarScan2 somatic analysis workflow |

Upon completion of the Somatic Variant Analysis, there is a `make somatic-completeion` command that is useful for collecting the run-related logs.
```
make somatic-completion
```

## Troubleshooting
If an error is encountered while deploying or using the pipeline, please open an [issue](https://github.com/pblaney/mgp1000/issues) so that it can be addressed and others who may have a similar issue can have a resource for potential solutions.

To further facilitate this, a cataloge of common issues and their solutions will be maintained within the [Wiki](https://github.com/pblaney/mgp1000/wiki)


## Contact
If there are any further questions, suggestions for improvement, or wishes for collaboration please feel free to email: patrick.blaney@nyulangone.org
