#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Script adapted and modified from Poelstra et al. 2021, Systematic Biology (https://doi.org/10.1093/sysbio/syaa053)

## Software:
# BEDtools needs to be included in $PATH (v2.30.0; https://bedtools.readthedocs.io/en/latest/)
# BEDOPS needs to be included in $PATH (v2.4.38; https://bedops.readthedocs.io/en/latest/)

## Command-line args:
vcf_altref=$1
vcf_filt_mask=$2
bed_removed_sites=$3
bed_dir=$4

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### 01_maskbed.sh: Starting script."
echo -e "#### 01_maskbed.sh: Raw VCF: $vcf_altref"
echo -e "#### 01_maskbed.sh: Fully filtered VCF: $vcf_filt_mask"
echo -e "#### 01_maskbed.sh: BED file for removed sites: $bed_removed_sites"
echo -e "#### 01_maskbed.sh: Directory for BED file: $bed_dir \n\n"

################################################################################
#### CREATE BED FILE WITH MASKED SITES ####
################################################################################
mkdir -p $bed_dir/tmpdir

echo -e "#### 01_maskbed.sh: Creating $bed_removed_sites with masked sites ..."
bedtools intersect -v -a <(vcf2bed --sort-tmpdir=$bed_dir/tmpdir < $vcf_altref | cut -f 1,2,3) -b <(vcf2bed --sort-tmpdir=$bed_dir/tmpdir < $vcf_filt_mask | cut -f 1,2,3) > $bed_removed_sites

## Report:
echo -e "\n#### 01_maskbed.sh: Line count of removed sites: $(wc -l < $bed_removed_sites)."
echo -e "#### 01_maskbed.sh: Done with script."
date

