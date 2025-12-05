#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################

## Command-line args:
scripts_dir=$1
in_file=$2
bam_hits=$3
set_id=$4

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### cov_plot.sh: Starting script."
echo -e "#### cov_plot.sh: Directory with scripts: $scripts_dir"
echo -e "#### cov_plot.sh: File with individuals: $in_file"
echo -e "#### cov_plot.sh: Directory with BAM hits: $bam_hits"
echo -e "#### cov_plot: Set ID: $set_id \n\n"

################################################################################
#### ESTIMATE AND PLOT COVERAGE DISTRIBUTIONS ####
################################################################################
echo -e "#### cov_plot.sh: Submitting R script to estimate and plot coverage distributions for each individual ...\n"
Rscript $scripts_dir/cov_plot.R $in_file $bam_hits $set_id

echo -e "\n#### cov_plot.sh: Done with script."
date






