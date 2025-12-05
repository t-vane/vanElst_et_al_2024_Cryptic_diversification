#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Script adapted and modified from Poelstra et al. 2021, Systematic Biology (https://doi.org/10.1093/sysbio/syaa053)

## Software:
# VCFtools needs to be included in $PATH (v0.1.17; https://vcftools.github.io/index.html)
# gatk needs to be included in $PATH (v4.1.9.0; https://gatk.broadinstitute.org/hc/en-us)

## Command-line args:
vcf_in=$1
vcf_out_soft=$2
vcf_out_hard=$3
reference=$4

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### 04_filter_gatk.sh: Starting script."
echo -e "#### 04_filter_gatk.sh: Input VCF: $vcf_in"
echo -e "#### 04_filter_gatk.sh: Output VCF for soft-filtering: $vcf_out_soft"
echo -e "#### 04_filter_gatk.sh: Output VCF for hard-filtering: $vcf_out_hard"
echo -e "#### 04_filter_gatk.sh: Reference genome: $reference \n\n"

################################################################################
#### SOFT-FILTER VCF FILE ####
################################################################################
echo -e "#### 04_filter_gatk.sh: Soft-filtering input VCF for\n"
echo -e "#### 04_filter_gatk.sh: FisherStrand, RMSMappingQuality, MappingQualityRankSumTest, ReadPosRankSumTest and AlleleBalance ...\n"
gatk VariantFiltration -R $reference -V $vcf_in -O $vcf_out_soft \
	--filter-expression "FS > 60.0" --filter-name "FS_gt60" \
	--filter-expression "MQ < 40.0" --filter-name "MQ_lt40" \
	--filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum_ltm12.5" \
	--filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum_ltm8" \
	--filter-expression "ABHet > 0.01 && ABHet < 0.2 || ABHet > 0.8 && ABHet < 0.99" --filter-name "ABhet_filt"
	
################################################################################
#### HARD-FILTER VCF FILE ####
################################################################################
echo -e "#### 04_filter_gatk.sh: Retaining only bi-allelic sites and hard-filtering input VCF for\n"
echo -e "#### 04_filter_gatk.sh: FisherStrand, RMSMappingQuality, MappingQualityRankSumTest, ReadPosRankSumTest and AlleleBalance ...\n"
vcftools --vcf $vcf_out_soft --remove-filtered-all --max-non-ref-af 0.99 --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --stdout > $vcf_out_hard

## Report:
nvar_in=$(grep -cv "^#" $vcf_in)
nvar_out=$(grep -cv "^#" $vcf_out_hard)
nvar_filt=$(( $nvar_in - $nvar_out ))

nfilt_fs=$(grep -c "FS_gt60" $vcf_out_soft || true)
nfilt_mq=$(grep -c "MQ_lt40" $vcf_out_soft || true)
nfilt_mqr=$(grep -c "MQRankSum_ltm12" $vcf_out_soft || true)
nfilt_readpos=$(grep -c "readposRankSum_ltm8" $vcf_out_soft || true)
nfilt_abhet=$(grep -c "ABhet_filt" $vcf_out_soft || true)

echo -e "\n\n"
echo -e "#### 04_filter_gatk.sh: Number of SNPs in input VCF: $nvar_in"
echo -e "#### 04_filter_gatk.sh: Number of SNPs in output VCF: $nvar_out"
echo -e "#### 04_filter_gatk.sh: Number of SNPs filtered: $nvar_filt\n"

echo -e "#### 04_filter_gatk.sh: Number of SNPs filtered by FS_gt60: $nfilt_fs"
echo -e "#### 04_filter_gatk.sh: Number of SNPs filtered by MQ_lt40: $nfilt_mq"
echo -e "#### 04_filter_gatk.sh: Number of SNPs filtered by MQRankSum_ltm12: $nfilt_mqr"
echo -e "#### 04_filter_gatk.sh: Number of SNPs filtered by readposRankSum_ltm8: $nfilt_readpos"
echo -e "#### 04_filter_gatk.sh: Number of SNPs filtered by ABhet_filt: $nfilt_abhet"

echo -e "\n#### 04_filter_gatk.sh: Listing output VCF:"
ls -lh $vcf_out_hard
[[ $(grep -cv "^#" $vcf_out_hard) = 0 ]] && echo -e "\n\n#### 04_filter_gatk.shh: ERROR: VCF is empty\n" >&2 && exit 1

echo -e "\n#### 04_filter_gatk.sh: Done with script."
date
