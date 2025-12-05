#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# NGSadmix needs to be included in $PATH (v32; http://www.popgen.dk/software/index.php/NgsAdmix)

## Command-line args:
nt=$1
k=$2
beagle=$3
out_dir=$4
minind=$5
set_id=$6
seed=$SLURM_ARRAY_TASK_ID

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### ngsadmix.sh: Starting script."
echo -e "#### ngsadmix.sh: Number of threads: $nt"
echo -e "#### ngsadmix.sh: Number of clusters: $k"
echo -e "#### ngsadmix.sh: Genotype likelihood file in beagle format: $beagle"
echo -e "#### ngsadmix.sh: Output directory: $out_dir"
echo -e "#### ngsadmix.sh: Minimum number of represented individuals: $minind"
echo -e "#### ngsadmix.sh: Set ID: $set_id"
echo -e "#### ngsadmix.sh: Seed (repetition number): $seed \n\n"

################################################################################
#### INFER INDIVIDUAL ANCESTRIES ####
################################################################################
echo -e "#### ngsadmix.sh: Estimating individual ancestries for $k clusters (repetition number $seed)...\n"
NGSadmix -P $nt -likes $beagle -seed $seed  -K $k -outfiles $out_dir/$set_id.K$k.seed$seed -minMaf 0.05 -minInd $minind -tol 0.000001

echo -e "\n#### ngsadmix.sh: Done with script."
date


