#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# Bcftools needs to be included in $PATH (v1.11; https://samtools.github.io/bcftools/)
# VCFtools needs to be included in $PATH (v0.1.17; https://vcftools.github.io/index.html)
gatk3=/home/nibtve93/software/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar # (v3.8.1; https://gatk.broadinstitute.org/hc/en-us)

## Command-line args:
reference=$1
vcf_in=$2
vcf_out=$3

## Activate conda environment
conda activate java

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### 01_filter_indels_invariants.sh: Starting script."
echo -e "#### 01_filter_indels_invariants.sh: Reference genome: $reference"
echo -e "#### 01_filter_indels_invariants.sh: Input VCF (compressed): $vcf_in"
echo -e "#### 01_filter_indels_invariants.sh: Output VCF: $vcf_out \n\n"

################################################################################
#### ANNOTATE WITH NS INFO FIELD AND REMOVE INDELS AND INVARIANT SITES ####
################################################################################
echo -e "#### 01_filter_indels_invariants.sh: Decompressing and annotating with number of samples with data (NS INFO field) ...\n"
bcftools +fill-tags $vcf_in -Ov -o $(dirname $vcf_in)/$(basename $vcf_in .vcf).annot.vcf -- -t NS

echo -e "#### 01_filter_indels_invariants.sh: Removing indels ...\n"
java -jar $gatk3 -T SelectVariants -R $reference -V $(dirname $vcf_in)/$(basename $vcf_in .vcf).annot.vcf -o $(dirname $vcf_in)/$(basename $vcf_in .vcf).tmp.vcf -selectType SNP
grep -v -E ',\*|\*,' $(dirname $vcf_in)/$(basename $vcf_in .vcf).tmp.vcf > $(dirname $vcf_in)/$(basename $vcf_in .vcf).tmp.tmp.vcf
mv $(dirname $vcf_in)/$(basename $vcf_in .vcf).tmp.tmp.vcf $(dirname $vcf_in)/$(basename $vcf_in .vcf).tmp.vcf

echo -e "#### 01_filter_indels_invariants.sh: Removing invariant sites ...\n"
vcftools --vcf $(dirname $vcf_in)/$(basename $vcf_in .vcf).tmp.vcf --recode --recode-INFO-all --max-non-ref-af 0.99 --min-alleles 2 --stdout > $vcf_out
rm $(dirname $vcf_in)/$(basename $vcf_in .vcf).tmp.vcf

## Report:
nvar_in=$(grep -cv "^#" $(dirname $vcf_in)/$(basename $vcf_in .vcf).annot.vcf)
nvar_out=$(grep -cv "^#" $vcf_out)
nvar_filt=$(( $nvar_in - $nvar_out ))

echo -e "\n\n"
echo -e "#### 01_filter_indels_invariants.sh: Number of sites in input VCF: $nvar_in"
echo -e "#### 01_filter_indels_invariants.sh: Number of sites in output VCF: $nvar_out"
echo -e "#### 01_filter_indels_invariants.sh: Number of sites filtered: $nvar_filt"
echo
echo -e "#### 01_filter_indels_invariants.sh: Listing output VCF:"
ls -lh $vcf_out
[[ $(grep -cv "^#" $vcf_out) = 0 ]] && echo -e "\n\n#### 01_filter_indels_invariants.sh: ERROR: VCF is empty\n" >&2 && exit 1

echo -e "\n#### 01_filter_indels_invariants.sh: Done with script."
date


