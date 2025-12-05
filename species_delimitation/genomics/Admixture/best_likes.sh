#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################

## Command-line args:
log=$1
like_file=$2
k=$3
seed=$4

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### print_likes.sh: Starting script."
echo -e "#### print_likes.sh: Respective log file: $log"
echo -e "#### print_likes.sh: Likelihoods summary file: $like_file"
echo -e "#### print_likes.sh: Number of clusters: $k"
echo -e "#### print_likes.sh: Seed: $seed \n\n"

################################################################################
#### EXTRACT LIKELIHOOD ####
################################################################################
echo -e "#### print_likes.sh: Extracting likelihood ...\n"
grep "best" $log | awk '{print $k}' | cut -d'=' -f2- | sort -g | sed "s/after/$k/g" | sed "s/iterations/$seed/g" >> $like_file

echo -e "#### print_likes.sh: Removing fopt.gz files ...\n"
rm $(dirname $log)/$(basename $log .log).fopt.gz

echo -e "\n#### print_likes.sh: Done with script."
date

