#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# gatk needs to be included in $PATH (v4.1.9.0; https://gatk.broadinstitute.org/hc/en-us)

## Command-line args:
vcf_list=$1
vcf_file=$2

## Activate conda environment
conda activate java

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### merge_vcfs.sh: Starting script."
echo -e "#### merge_vcfs.sh: List of per-scaffold VCF files: $vcf_list"
echo -e "#### merge_vcfs.sh: Output VCF file: $vcf_file \n\n"

################################################################################
#### MERGE PER-SCAFFOLD VCF FILES ####
################################################################################
echo -e "#### merge_vcfs.sh: Merging per-scaffold VCF files ...\n"
gatk MergeVcfs -I $vcf_list -O $vcf_file

echo -e "#### merge_vcfs.sh: Extracting VCF file ...\n"
gzip -d $vcf_file

## Report:
echo -e "\n#### merge_vcfs.sh: Done with script."
date







