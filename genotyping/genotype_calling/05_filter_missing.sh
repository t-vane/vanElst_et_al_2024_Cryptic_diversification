#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Script adapted and modified from Poelstra et al. 2021, Systematic Biology (https://doi.org/10.1093/sysbio/syaa053)

## Software:
# VCFtools needs to be included in $PATH (v0.1.17; https://vcftools.github.io/index.html)

## Command-line args:
vcf_in=$1
vcf_out=$2
maxmiss_geno=$3
filter_inds=$4
maxmiss_ind=$5

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### 05_filter_missing.sh: Starting script."
echo -e "#### 05_filter_missing.sh: Input VCF: $vcf_in"
echo -e "#### 05_filter_missing.sh: Output VCF: $vcf_out"
echo -e "#### 05_filter_missing.sh: Maximum missingness value for genotypes: $maxmiss_geno"
echo -e "#### 05_filter_missing.sh: NOTE: Maximum missingness is inverted for genotypes, with a value of 1 indicating no missing data (this is how the '--max-missing' option in VCFtools works)."
echo -e "#### 05_filter_missing.sh: Filter individuals: $filter_inds"
[[ $filter_inds ]] && echo -e "#### 05_filter_missing.sh: Maximum missingness value for individuals: $maxmiss_ind \n\n"

################################################################################
#### FILTER FOR MISSINGNESS PER GENOTYPE (AND INDIVIDUAL) ####
################################################################################
## Estimate number of indidviduals in input
echo -e "#### 05_filter_missing.sh: Estimating number of individuals in input VCF ...\n"
nind_in=$(bcftools query -l $vcf_in | wc -l || true)

## Create temporary filter files
mkdir -p $(dirname $vcf_out)/tmp
filterfile_prefix=$(dirname $vcf_out)/tmp/filterfile.$RANDOM
vcf_tmp=$(dirname $vcf_out)/tmp/vcf_indfilter.$RANDOM.vcf

## Run filtering
if [[ $filter_inds ]]
then
	echo -e "#### 05_filter_missing.sh: Filtering individuals by missing data ... "
	# Get amount of missing data per individual
	vcftools --vcf $vcf_in --missing-indv --stdout > $filterfile_prefix.imiss
	# Get list of individuals with too much missing data
	tail -n +2 $filterfile_prefix.imiss | awk -v var=$maxmiss_ind '$5 > var' | cut -f1 > $filterfile_prefix.HiMissInds
	# Remove individuals with too much missing data
	vcftools --vcf $vcf_in --remove $filterfile_prefix.HiMissInds --recode --recode-INFO-all --stdout > $vcf_tmp	
else
	echo -e "#### 05_filter_missing.sh: Only filtering by missing data at genotype level (no individuals will be removed)"
	vcf_tmp=$vcf_in
fi

echo -e "#### 05_filter_missing.sh: Filtering genotypes by missing data... "
vcftools --vcf $vcf_tmp --max-non-ref-af 0.99 --min-alleles 2 --max-missing $maxmiss_geno --recode --recode-INFO-all --stdout > $vcf_out

## Remove temporary files
rm $filterfile_prefix*
rm $vcf_tmp

## Report:
nvar_in=$(grep -cv "^#" $vcf_in || true)
nvar_out=$(grep -cv "^#" $vcf_out || true)
nvar_filt=$(( $nvar_in - $nvar_out ))

nind_filt=$(wc -l < $filterfile_prefix.HiMissInds || true)
nind_out=$(bcftools query -l $vcf_out | wc -l || true)

if [[ $filter_inds ]]
then
	echo -e "\n#### 05_filter_missing.sh: Number of indidviduals before individual filtering: $nind_in"
	echo -e "#### 05_filter_missing.sh: Number of indidviduals filtered: $nind_filt"
	echo -e "#### 05_filter_missing.sh: Number of indidviduals after individual filtering: $nind_out"
fi

echo -e "\n#### 05_filter_missing.sh: Number of SNPs before genotype filtering: $nvar_in"
echo -e "#### 05_filter_missing.sh: Number of SNPs filtered: $nvar_filt"
echo -e "#### 05_filter_missing.sh: Number of SNPs after genotype filtering: $nvar_out"

echo -e "\n#### 05_filter_missing.sh: Listing output VCF:"
ls -lh $vcf_out
[[ $(grep -cv "^#" $vcf_out) = 0 ]] && echo -e "\n\n#### 05_filter_missing.sh: ERROR: VCF is empty\n" >&2 && exit 1

echo -e "\n#### 05_filter_missing.sh: Done with script."
date

