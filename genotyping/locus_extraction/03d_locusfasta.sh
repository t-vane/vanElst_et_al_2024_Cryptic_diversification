#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Script adapted and modified from Poelstra et al. 2021, Systematic Biology (https://doi.org/10.1093/sysbio/syaa053)

## Software:
# SAMtools needs to be included in $PATH (v1.11; http://www.htslib.org/)

## Command-line args:
locuslist=$1
locusfasta_dir_intermed=$2
fasta_merged=$3
locus=$(sed -n "$SLURM_ARRAY_TASK_ID"p $locuslist)

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### 03d_locusfasta.sh: Starting script."
echo -e "#### 03d_locusfasta.sh: List of loci: $locuslist"
echo -e "#### 03d_locusfasta.sh: Directory for intermediate locus FASTA files: $locusfasta_dir_intermed"
echo -e "#### 03d_locusfasta.sh: Merged FASTA file: $fasta_merged"
echo -e "#### 03d_locusfasta.sh: Locus: $locus \n\n"

################################################################################
#### CREATE BY-LOCUS FASTA FILE  ####
################################################################################
echo -e "\n#### 03d_locusfasta.sh: Creating FASTA file for locus $locus ..."
locus_faidx=$(echo $locus | sed 's/:/,/g')
fasta=$locusfasta_dir_intermed/$locus.fa
faidx --regex $locus_faidx $fasta_merged | sed 's/,/:/g' > $fasta

## Report:
echo -e "\n#### 03d_locusfasta.sh: Done with script."
date

