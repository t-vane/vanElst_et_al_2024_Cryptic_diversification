#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# Bcftools needs to be included in $PATH (v1.11; https://samtools.github.io/bcftools/)
# VCFtools needs to be included in $PATH (v0.1.17; https://vcftools.github.io/index.html)

## Command-line args:
mindepthind=$1
maxdepthind=$2
gmin=$3
gmax=$4
minind=$5
mac=$6
vcf_in=$7
vcf_out=$8

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### 02_filter_dp.sh: Starting script."
echo -e "#### 02_filter_dp.sh: Minimum depth per individual: $mindepthind"
echo -e "#### 02_filter_dp.sh: Maximum depth per individual: $maxdepthind"
echo -e "#### 02_filter_dp.sh: Minimum global depth: $gmin"
echo -e "#### 02_filter_dp.sh: Maximum global depth: $gmax"
echo -e "#### 02_filter_dp.sh: Minimum number of individuals: $minind"
echo -e "#### 02_filter_dp.sh: Minor allele count: $mac"
echo -e "#### 02_filter_dp.sh: Input VCF: $vcf_in"
echo -e "#### 02_filter_dp.sh: Output VCF: $vcf_out \n\n"

################################################################################
#### FILTER BASED ON DEPTH, MINIMUM NUMBER OF INDIVIDUALS AND MINOR ALLELE COUNT ####
################################################################################
echo -e "#### 02_filter_dp.sh: Filtering based on depth and minimum number of individuals ...\n"
bcftools filter -e "INFO/DP<$gmin || INFO/DP>$gmax || INFO/NS<$minind" $vcf_in | bcftools filter -S . -e  "FORMAT/DP<$mindepthind" | bcftools filter -S . -e  "FORMAT/DP>$maxdepthind" -Ov -o $(dirname $vcf_in)/$(basename $vcf_in .vcf).tmp.vcf

echo -e "#### 02_filter_dp.sh: Filtering based on minor allele count ...\n"
vcftools --vcf $(dirname $vcf_in)/$(basename $vcf_in .vcf).tmp.vcf --mac $mac --recode --recode-INFO-all --stdout > $vcf_out
rm $(dirname $vcf_in)/$(basename $vcf_in .vcf).tmp.vcf

## Report:
nvar_in=$(grep -cv "^#" $vcf_in)
nvar_out=$(grep -cv "^#" $vcf_out)
nvar_filt=$(( $nvar_in - $nvar_out ))

echo -e "\n\n"
echo -e "#### 02_filter_dp.sh: Number of sites in input VCF: $nvar_in"
echo -e "#### 02_filter_dp.sh: Number of sites in output VCF: $nvar_out"
echo -e "#### 02_filter_dp.sh: Number of sites filtered: $nvar_filt"
echo
echo -e "#### 02_filter_dp.sh: Listing output VCF:"
ls -lh $vcf_out
[[ $(grep -cv "^#" $vcf_out) = 0 ]] && echo -e "\n\n#### 02_filter_dp.sh: ERROR: VCF is empty\n" >&2 && exit 1

echo -e "\n#### 02_filter_dp.sh: Done with script."
date


