#!/usr/bin/env bash
# This script will generate a SLURM batch submission script for any pipeline run command and
# then submit this to the scheduler

####################	Help Message	####################
Help()
{
	# Display help message
	echo "This script will generate a SLURM batch submission script for any pipeline"
	echo "run command and then submit this to the scheduler"
	echo ""
	echo "Syntax:"
	echo '	./nextflow_run_slurm_submitter.sh [-h] [pipelineStepScript] [runId] [userEmail]'
	echo '	[estimatedRuntime] "[moduleLoadCmd]" "[pipelineStepOptions]" "[additionalSlurmOptions]"'
	echo ""
	echo "Argument Descriptions:"
	echo "	[-h]				Print this message"
	echo "	[pipelineStepScript]		The name of the Nextflow script that will be run"
	echo "	[runId]				The unique ID associated with the run"
	echo "	[userEmail]			The email address to be sent notifications of pipeline progress"
	echo "	[estimatedRuntime]		The estimated number of days for run to complete"
	echo "	[moduleLoadCmd]			The text string that will be used to load required modules"
	echo "                 			within the parent SLURM job"
	echo "	[pipelineStepOptions]		The text string of any pipeline step specific run options"
	echo "	[additionalSlurmOptions]	OPTIONAL; The text string of any unqiue additional SLURM"
	echo "	                        	submission options"
	echo ""
}

while getopts ":h" option;
	do
		case $option in
			h) # Show help message
				Help
				exit;;
		   \?) # Reject other passed options
				echo "Invalid option"
				exit;;
		esac
	done

############################################################

# Capture command line arguments
pipelineStepScript=$1
runId=$2
userEmail=$3
estimatedRuntime=$4
moduleLoadCmd=$5
pipelineStepOptions=$6
additionalSlurmOptions=${7:-""}

# Create SLURM batch submission script with user-defined parameters
pipelineStep=$(echo "${pipelineStepScript}" | sed 's|.nf||')
submissionScript="slurmsub.${pipelineStep}.${runId}.sh"

echo "#!/bin/bash" > "${submissionScript}"
echo "" >> "${submissionScript}"
echo "#SBATCH --mail-user=${userEmail}" >> "${submissionScript}"
echo "#SBATCH --mail-type=BEGIN,END,FAIL" >> "${submissionScript}"
echo "#SBATCH --nodes=1" >> "${submissionScript}"
echo "#SBATCH --ntasks-per-node=1" >> "${submissionScript}"
echo "#SBATCH --cpus-per-task=2" >> "${submissionScript}"
echo "#SBATCH --mem=4G" >> "${submissionScript}"
echo "#SBATCH --time=${estimatedRuntime}-00:00:00" >> "${submissionScript}"
echo "#SBATCH --export=ALL" >> "${submissionScript}"
echo "#SBATCH --job-name=slurmsub.${pipelineStep}.${runId}" >> "${submissionScript}"
echo "#SBATCH --output=slurmsub.${pipelineStep}.${runId}.txt" >> "${submissionScript}"
echo "#SBATCH --error=slurmsub.${pipelineStep}.${runId}.err" >> "${submissionScript}"

if [[ "${additionalSlurmOptions}" != "" ]]; then
	echo "#SBATCH ${additionalSlurmOptions}" >> "${submissionScript}"
fi

echo "" >> "${submissionScript}"
echo "### ^^^   SLURM sbatch options   ^^^ ###" >> "${submissionScript}"
echo "############################################################" >> "${submissionScript}"
echo "" >> "${submissionScript}"

echo "eval ${moduleLoadCmd}" >> "${submissionScript}"

echo "" >> "${submissionScript}"
echo "### ^^^   Modules loaded   ^^^ ###" >> "${submissionScript}"
echo "############################################################" >> "${submissionScript}"
echo "" >> "${submissionScript}"

echo "export NXF_ANSI_LOG=false" >> "${submissionScript}"
echo 'export NXF_OPTS="-Xms500M -Xmx2G -Dleveldb.mmap=false"' >> "${submissionScript}"

runCommand="
./nextflow run \
-resume \
${pipelineStepScript} \
--run_id ${runId} \
--email ${userEmail} \
-profile ${pipelineStep} \
${pipelineStepOptions}
"
echo "${runCommand}" >> "${submissionScript}"

echo "### ^^^   Nexflow run options and command   ^^^ ###" >> "${submissionScript}"
echo "############################################################" >> "${submissionScript}"

chmod +x "${submissionScript}"

eval sbatch ${submissionScript}
