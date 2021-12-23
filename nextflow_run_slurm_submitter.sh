#!/bin/bash
# This script will generate a SLURM batch submission script for any pipeline run command and
# then submit this to the scheduler

# Capture command line arguments
pipelineStepScript=$1
runId=$2
userEmail=$3
estimatedRuntime=$4
moduleLoadCmd=$5
pipelineStepOptions=$6

# Create SLURM batch submission script with user-defined parameters
pipelineStep=$(echo "${pipelineStepScript}" | sed 's|.nf||')
submissionScript="slurmsub.${pipelineStep}.${runId}.sh"

runtime=$((${estimatedRuntime} * 24))

touch "${submissionScript}"
chmod +x "${submissionScript}"

echo "#!/bin/bash" >> "${submissionScript}"
echo "" >> "${submissionScript}"
echo "#SBATCH --mail-user=${userEmail}" >> "${submissionScript}"
echo "#SBATCH --mail-type=BEGIN,END,FAIL" >> "${submissionScript}"
echo "#SBATCH --ntasks-per-node=1" >> "${submissionScript}"
echo "#SBATCH --cpus-per-task=2" >> "${submissionScript}"
echo "#SBATCH --mem=4G" >> "${submissionScript}"
echo "#SBATCH --time=${runtime}:00:00" >> "${submissionScript}"
echo "#SBATCH --job-name=slurmsub.${pipelineStep}.${runId}" >> "${submissionScript}"
echo "#SBATCH --output=slurmsub.${pipelineStep}.${runId}.txt" >> "${submissionScript}"
echo "#SBATCH --error=slurmsub.${pipelineStep}.${runId}.err" >> "${submissionScript}"

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
echo 'export NXF_OPTS="-Xms500M -Xmx2G -Dleveldb.mmap=false"'
echo ""

runCommand="
nextflow run \
${pipelineStepScript} \
--run_id ${runId} \
--email ${userEmail} \
-profile ${pipelineStep} \
${pipelineStepOptions}
"
echo "${runCommand}" >> "${submissionScript}"

echo "### ^^^   Nexflow run options and command   ^^^ ###" >> "${submissionScript}"
echo "############################################################" >> "${submissionScript}"

eval sbatch ${submissionScript}
