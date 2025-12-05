#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Script adapted and modified from Poelstra et al. 2021, Systematic Biology (https://doi.org/10.1093/sysbio/syaa053)

## Software:
# BEDtools needs to be included in $PATH (v2.30.0; https://bedtools.readthedocs.io/en/latest/)

## Command-line args:
locusbed_intermed=$1
locusbed_final=$2
vcf_highdepth=$3
vcf_filt_intersect=$4

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### 03b_intersect.sh: Starting script."
echo -e "#### 03b_intersect.sh: Intermediate locus BED file: $locusbed_intermed"
echo -e "#### 03b_intersect.sh: Final locus BED file: $locusbed_final"
echo -e "#### 03b_intersect.sh: VCF file with loci with too high depth: $vcf_highdepth"
echo -e "#### 03b_intersect.sh: Fully filtered VCF file: $vcf_filt_intersect \n\n"

################################################################################
#### INTERSECT BED FILE WITH LOCI WITH TOO HIGH DEPTH ####
################################################################################
echo -e "#### 03b_intersect.sh: Intersecting BED file with loci with too high depth ..."
bedtools intersect -v -a $locusbed_intermed -b $vcf_highdepth > $locusbed_final

## Report:
echo -e "#### 03b_intersect.sh: Number of loci after removing SNPs with too high depth: $(wc -l < $locusbed_final)"

snps_in_loci=$(bedtools intersect -u -a $vcf_filt_intersect -b $locusbed_final | grep -cv "##")
snps_in_vcf=$(grep -cv "##" $vcf_filt_intersect)
snps_lost=$(( $snps_in_vcf - $snps_in_loci ))

echo -e "#### 03b_intersect.sh: Number of SNPs in VCF: $snps_in_vcf"
echo -e "#### 03b_intersect.sh: Number of SNPs in loci: $snps_in_loci"
echo -e "#### 03b_intersect.sh: Number of lost SNPs (in VCF but not in loci): $snps_lost"

echo -e "\n#### 03b_intersect.sh: Done with script."
date

