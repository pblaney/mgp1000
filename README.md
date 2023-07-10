<div align="center">
	<img alt="MGP1000 logo" src="https://raw.githubusercontent.com/pblaney/mgp1000/4c6c86d956b30e6b64bdad50619c9f5b76cee2d9/docs/mgp1000Logo.svg" width="250" />

# Myeloma Genome Pipeline 1000
Comprehensive bioinformatics pipeline for the large-scale collaborative analysis of Multiple Myeloma genomes in an effort to delineate the broad spectrum of somatic events
</div>

## Pipeline Overview
In order to analyze over one thousand matched tumor/normal whole-genome samples across multiple data centers in a consistent manner, a pipeline was created that leverages the workflow management, portability, and reproducibility of [Nextflow](http://www.nextflow.io/) in conjuction with [Singularity](https://sylabs.io/docs/).

The entire pipeline is divided into 3 modules: Preprocessing, Germline, and Somatic
<div align="center">
	<img alt="Pipeline flowchart" src="https://github.com/pblaney/mgp1000/blob/development/docs/pipelineArchitectureForGitHub.png?raw=true" width="475" />
</div>

## Deploying the Pipeline
The pipeline was developed to be run on various HPCs without concern of environment incompatibilities, version issues, or missing dependencies. None of the commands require admin access or `sudo` to be completed. However, there are a few assumptions regarding initial setup of the pipeline but the required software should be readily available on nearly all HPC systems.
* Git
* GNU Utilities
* Java 11 (or later)
* Singularity (validated on v3.1, v3.5.2, v3.7.1, v3.9.8 other versions will be tested)


## Clone GitHub Repository
The first step in the deployment process is to clone the MGP1000 GitHub repository to a location on your HPC that is large enough to hold the input/output data and has access to the job scheduling software, such as SLURM or SGE.
```
git clone https://github.com/pblaney/mgp1000.git
```

### Installing Git LFS
In an effort to containerize the pipeline further, all the necessary reference files and Singularity container images are stored in the GitHub repository using their complementary [Large File Storage (LFS)](https://git-lfs.github.com) extension. This requires a simple installation of the binary executible file at a location on your `$PATH`. The extension pairs seemlessly with Git to download all files while cloning the repository.

**NOTE: Many HPC environments may already have this dependency installed, if so users only need to run the last two code blocks.**

If required, a `make` command will complete the installation of Linux AMD64 binary executible git-lfs file (v3.2.0). Other binary files available [here](https://github.com/git-lfs/git-lfs/releases)


```
make install-gitlfs-linuxamd64
```
Move the `git-lfs` binary to a location on `$PATH`
```
mv git-lfs $HOME/bin
```
Complete the install and configuration
```
git-lfs install
```
Use `git-lfs` to complete the clone
```
git-lfs pull
```

### Reference Data
To facilitate ease of use, reproducibility, and consistency between all users of the pipeline, all required reference data has been provided within the `references/hg38/` directory. Detailed provenance of each file per tool is included in the pipline [Wiki](https://github.com/pblaney/mgp1000/wiki) for full traceability.

### Containers
For the same reasons as with the reference data, the Singularity image files needed for each tool's container is provided within the `containers/` directory. All containers were originally developed with Docker and all tags can be found on the associated [DockerHub](https://hub.docker.com/r/patrickblaneynyu/mgp1000).

## Install Nextflow
This series of commands will install Nextflow. The most current version of Nextflow requires Java 11 or later. Therefore, the user may need to load this version of Java to complete the install.

**NOTE: Many HPC environments may already have this dependency installed, if so this section can be skipped.**

```
java -version
```
```
make install-nextflow
```
For ease of use, ove the binary executible `nextflow` file to same directory as `git-lfs`
```
mv nextflow $HOME/bin
```


## Prepare the Pipeline for Usage
Due to size, certain reference files are GNU zipped so the `make prep-pipeline` command must be run to prepare them for use in the pipeline. Additionally, an `input` directory is created for staging all input BAM or FASTQ files, `preprocessedBams` subdirectory for BAMs that have undergone preprocessing and are ready for Germline/Somatic modules, and a `logs` directory to store Nextflow output log files for each run.
```
make prep-pipeline
```


## Stage Input BAM or FASTQ Files
By default, all input files are handled out of the `input` and `input/preprocessedBams` directories for the Preprocessing and Germline/Somatic modules, respectively. However, each module in the pipeline includes an option (`--input_dir`) for the user to define the input directory. Additionally, the pipeline will follow symbolic links for input files so there is no need to move files for staging. Given the possible size of the input data, the samples may have to be processed in batches. Additionally, the pipeline is designed to process batches of identical format, i.e. all BAMs or all FASTQs.

The Preprocessing module also supports lane split FASTQs as input. The files will be merged internally as part of the module run and will not alter any input files.

**NOTE: A key assumption is that any input FASTQs use an 'R1/R2' naming convention to designate paired-end read files.**

Example of staging input data files with symbolic link
```
ln -s /absolute/path/to/unprocessed/samples/directory/*.fastq.gz input/
```


## Run the Preprocessing Module
The Preprocessing module of the pipeline will be started with one command that will handle linking each individual process in the pipeline to the next. A key advantage of using Nextflow within an HPC environment is that will also perform all the job scheduling/submitting given the correct configuration with the user's [executor](https://www.nextflow.io/docs/latest/executor.html).

Currently supported executors:
* SLURM
* LSF

There are two methods for running each module in the pipeline: directly from the current environment or batch submission. An example of direct submission is given below and an example of batch submission is provided within the [Wiki](https://github.com/pblaney/mgp1000/wiki/Usage).

### Direct Submission Example
First, the user will need to load the necessary software to run the pipeline module to the environment. At most, this will require Java, Nextflow, and Singularity. This is a user-specific step so the commands may be different depending on the user's HPC configuration.

```
nextflow run preprocessing.nf \
  -bg \
  --run_id batch1 \
  --input_format fastq \
  --email user@example.com \
  -profile preprocessing
```

Here is the full help message for the Preprocessing module.
```
$ nextflow run preprocessing.nf --help

Usage:
  nextflow run preprocessing.nf --run_id STR --input_format STR -profile preprocessing
  [-bg] [-resume] [--lane_split STR] [--input_dir PATH] [--output_dir PATH] [--email STR]
  [--seq_protocol STR][--cpus INT] [--memory STR] [--queue_size INT] [--executor STR] [--help]

Mandatory Arguments:
  --run_id                       STR  Unique identifier for pipeline run
  --input_format                 STR  Format of input files
                                      [Default: fastq | Available: fastq, bam]
  -profile                       STR  Configuration profile to use, must use preprocessing                               

Main Options:
  -bg                           FLAG  Runs the pipeline processes in the background, this
                                      option should be included if deploying pipeline with
                                      real data set so processes will not be cut if user
                                      disconnects from deployment environment
  -resume                       FLAG  Successfully completed tasks are cached so that if
                                      the pipeline stops prematurely the previously
                                      completed tasks are skipped while maintaining their
                                      output
  --lane_split                   STR  Determines if input FASTQs are lane split per R1/R2
                                      [Default: no | Available: yes, no]
  --input_dir                   PATH  Directory that holds BAMs and associated index files,
                                      this should be given as an absolute path
                                      [Default: input/]
  --output_dir                  PATH  Directory that will hold all output files this should
                                      be given as an absolute path
                                      [Default: output/]
  --email                        STR  Email address to send workflow completion/stoppage
                                      notification
  --seq_protocol                 STR  Sequencing protocol of the input, WGS for whole-genome
                                      and WES for whole-exome
                                      [Default: WGS | Available: WGS, WES]
  --cpus                         INT  Globally set the number of cpus to be allocated
  --memory                       STR  Globally set the amount of memory to be allocated,
                                      written as '##.GB' or '##.MB'
  --queue_size                   INT  Set max number of tasks the pipeline will launch
                                      [Default: 100]
  --executor                     STR  Set the job executor for the run
                                      [Default: slurm | Available: local, slurm, lsf]
  --help                        FLAG  Prints this message
```


### Preprocessing Output Description
This module of the pipeline generates various per tool QC metrics that are useful in determining samples for best use in downstream analyses. By default, all output files are stored into process-specific subdirectories within the `output/` directory.

Here is a snapshot of the expected subdirectories within the `output/preprocessing` base directory after a successful run of the [Preprocessing module](https://github.com/pblaney/mgp1000/wiki/Preprocessing):

| Subdirectory | Output Files | Description of Files |
| --- | --- | --- |
| `trimLogs` | `*.trim.log` | number of reads before and after trimming for quality |
| `fastqc` | `*_fastqc.[html / zip]` | in-depth quality evaluation on a per base and per sequence manner | 
| `alignmentFlagstats` | `*.alignment.flagstat.log` | initial number of reads of various alignment designation |
| `markdupFlagstats` | `*.markdup.[log / flagstat.log]` | number of detected duplicate reads, number of reads after deduplication |
| `finalPreprocessedBams` | `*.final.[bam / bai]` | final preprocessed BAM and index for downstream analysis |
| `coverageMetrics` | `*.coverage.metrics.txt` | genome-wide coverage metrics of final BAM |
| `gcBiasMetrics` | `*.gcbias.[metrics.txt / metrics.pdf / summary.txt]` | genome-wide GC bias metrics of final BAM |

Upon completion of the Preprocessing run, there is a `make preprocessing-completion` command that is useful for collecting the run-related logs.
```
make preprocessing-completion
```


## Run the Germline Module
The most important component of this module of the pipeline is the user-provided sample sheet CSV. This file includes two comma-separated columns: filename of normal sample BAMs and filename of corresponding paired tumor sample BAMs. An example of this is provided in `samplesheet.csv` within the `testSamples` directory. The sample sheet file should typically be within the main `mgp1000` directory.

```
nextflow run germline.nf \
  -bg \
  --run_id batch1 \
  --sample_sheet wgs_samples.csv \
  --email user@example.com \
  -profile germline
```

Here is the full help message for the Germline module.
```
$ nextflow run germline.nf --help

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
```

### Germline Output Description
This module of the pipeline generates a per-patient VCF and admixture estimation of individual ancestries in the context of 23 populations outlined in the 1000 Genomes Project. By default, all output files are stored within the `output/` directory.

Here is a snapshot of the expected subdirectories within the `output/germline/` base directory after a successful run of the [Germline module](https://github.com/pblaney/mgp1000/wiki/Germline):

| Germline Output File | Description of File |
| --- | --- |
| `*.deepvariant.germline.snp.indel.vcf.gz` | germline SNP and InDel VCF |
| `*.fastngsadmix.23.qopt` | per patient estimated admixture proportions, first two rows with the names of the populations analyzed and the converged upon estimates, and then 100 rows with the bootstrapping estimates |

Upon completion of the Germline module there is a `make germline-completion` command that is useful for collecting the run-related logs.
```
make germline-completion
```


## Run the Somatic Module
This module uses the same user-provided sample sheet CSV as the Germline module.

### Note on Parameters
There is one parameter that will prepare necessary reference file as part of this module of the pipeline: `--mutect_ref_vcf_concatenated`. This parameters only needs be set to `no` for the **first run** of the Somatic module of the pipeline. By default, all tools in this module will be used with the standard command in the usage example. The consensus output per variant type expects all tools to be included in the run for consensus processes to be run.
```
nextflow run somatic.nf --help

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
```

### Somatic Output Description
This module of the pipeline generates per tumor-normal pair consensus calls for SNVs, InDels, CNVs, and SVs, with some additional supporting tools for telomere analysis and graph-based inference of rearrangement state. Each tool used has its native output kept within a self-named subdirectory while the final consensus output files per tumor-normal pair are funneled into the `consensus` subdirectory.

Here is a snapshot of the final `output/somatic/consensus/[tumor_normal_id]` directory after a successful run of the [Somatic module](https://github.com/pblaney/mgp1000/wiki/Somatic):

| Somatic Consensus Output File | Description of File |
| --- | --- |
| `*.hq.union.consensus.somatic.snv.txt.gz` | consensus SNV table |
| `*.hq.union.consensus.somatic.indel.txt.gz` | consensus InDel table |
| `*.hq.union.consensus.somatic.cnv.bed` | consensus CNV BED |
| `*.hq.union.consensus.somatic.sv.bedpe` | consensus SV BEDPE |
| `*.hq.union.consensus.somatic.sv.intersection.plot.pdf` | intersection of SV caller consensus |

Additional per-tool subdirectories included in the base `output/somatic` output directory:

| Subdirectory | Description of Files |
| --- | --- |
| `battenberg` | native output of Battenberg |
| `conpair` | native output of Conpair |
| `delly` | native output of DELLY2 |
| `facets` | native output of FACETS |
| `igcaller` | native output of IgCaller |
| `manta` | native output of Manta |
| `mutect` | native output of Mutect2 |
| `sexOfSamples` | sample sex estimation using alleleCount |
| `strelka` | native output of Strelka2 |
| `svaba` | native output of SvABA |
| `telomerecat` | native output of Telomerecat |
| `telomereHunter` | native output of TelomereHunter |
| `varscan` | native output of VarScan2 |

Upon completion of the Somatic module, there is a `make somatic-completion` command that is useful for collecting the run-related logs.
```
make somatic-completion
```

## Troubleshooting
If an error is encountered while deploying or using the pipeline, please open an [issue](https://github.com/pblaney/mgp1000/issues) so that it can be addressed and others who may have a similar issue can have a resource for potential solutions.

To further facilitate this, a cataloge of common issues and their solutions will be maintained within the [Wiki](https://github.com/pblaney/mgp1000/wiki)


## Citation
Hopefully soon....


## Contact
If there are any further questions, suggestions for improvement, or wishes for collaboration please feel free to email: patrick.blaney@nyulangone.org
