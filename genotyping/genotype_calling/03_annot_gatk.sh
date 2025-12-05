#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Script adapted and modified from Poelstra et al. 2021, Systematic Biology (https://doi.org/10.1093/sysbio/syaa053)

## Software:
# SAMtools needs to be included in $PATH (v1.11; http://www.htslib.org/)
# VCFtools needs to be included in $PATH (v0.1.17; https://vcftools.github.io/index.html)
gatk3=/home/nibtve93/software/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar # (v3.8.1; https://gatk.broadinstitute.org/hc/en-us)

## Command-line args:
vcf_in=$1
vcf_out=$2
bam_dir=$3
reference=$4
suffix=$5

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### 03_annot_gatk.sh: Starting script."
echo -e "#### 03_annot_gatk.sh: Input VCF: $vcf_in"
echo -e "#### 03_annot_gatk: Output VCF: $vcf_out"
echo -e "#### 03_annot_gatk: Directory with BAM files: $bam_dir"
echo -e "#### 03_annot_gatk.sh: Reference genome: $reference"
echo -e "#### 03_annot_gatk.sh: Suffix for BAM files: $suffix \n\n"

################################################################################
#### ANNOTATE VCF FILE ####
################################################################################
## Sort VCF file
echo -e "#### 03_annot_gatk.sh: Sorting $vcf_in ...\n"
vcf-sort $vcf_in > $vcf_in.tmp
mv $vcf_in.tmp > $vcf_in

## Create list of BAM files
echo -e "#### 03_annot_gatk.sh: Creating list of BAM files for individuals in $vcf_in ...\n"
for i in $(bcftools query -l $vcf_in)
do
	echo "$bam_dir/${i}.$suffix.bam"
done > $(dirname $vcf_in)/bamFiles.txt

## Annotate VCF file
# Process list of BAM files to create command argument
bam_arg=""
while read -r bam
do
  [[ ! -f $bam.bai ]] && samtools index $bam
  bam_arg=$(echo "$bam_arg -I $bam")
done < $(dirname $vcf_in)/bamFiles.txt

# Run GATK for annotation
echo -e "#### 03_annot_gatk.sh: Annotating vcf_in with INFO fields for\n"
echo -e "#### 03_annot_gatk.sh: FisherStrand, RMSMappingQuality, MappingQualityRankSumTest, ReadPosRankSumTest and AlleleBalance ...\n"
java -jar $gatk3 -T VariantAnnotator -R $reference -V $vcf_in -o $vcf_out $bam_arg -A FisherStrand -A RMSMappingQuality -A MappingQualityRankSumTest -A ReadPosRankSumTest -A AlleleBalance

## Report:
echo -e "\n#### 03_annot_gatk.sh: Listing output VCF:"
ls -lh $vcf_out
[[ $(grep -cv "^#" $vcf_out) = 0 ]] && echo -e "\n\n#### 03_annot_gatk.sh: ERROR: VCF is empty\n" >&2 && exit 1

echo -e "\n#### 03_annot_gatk.sh: Done with script."
date

