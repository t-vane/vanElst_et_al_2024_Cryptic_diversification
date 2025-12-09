#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# SAMtools needs to be included in $PATH (v1.11; http://www.htslib.org/)
# gatk needs to be included in $PATH (v4.1.9.0; https://gatk.broadinstitute.org/hc/en-us)

## Command-line args:
nt=$1
mem=$2
reference=$3
ind_file=$4
bam_dir=$5
suffix=$6
gvcf_dir=$7

indv=$(sed -n "$SLURM_ARRAY_TASK_ID"p $ind_file)

## Activate conda environment
conda activate java

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### haplotype_caller.sh: Starting script."
echo -e "#### haplotype_caller.sh: Number of threads: $nt"
echo -e "#### haplotype_caller.sh: Memory: $mem"
echo -e "#### haplotype_caller.sh: Reference genome: $reference"
echo -e "#### haplotype_caller.sh: List with individuals: $ind_file"
echo -e "#### haplotype_caller.sh: Directory with BAM files: $bam_dir"
echo -e "#### haplotype_caller.sh: Suffix of BAM files: $suffic"
echo -e "#### haplotype_caller.sh: Output directory for GVCF files: $gvcf_dir"
echo -e "#### haplotype_caller.sh: Individual: $indv \n\n"

################################################################################
#### CREATE GVCF FILE ####
################################################################################
echo -e "#### haplotype_caller.sh: Indexing BAM file for individual $indv ...\n"
[[ ! -f $bam_dir/$indv.$suffix.bai ]] && samtools index $bam_dir/$indv.$suffix

echo -e "#### haplotype_caller.sh: Creating GVCF file for individual $indv ...\n"
gatk --java-options "-Xmx${mem}g" HaplotypeCaller -R $reference -I $bam_dir/$indv.$suffix -O $gvcf_dir/$indv.rawvariants.g.vcf.gz -ERC GVCF --pairHMM AVX_LOGLESS_CACHING_OMP --native-pair-hmm-threads $nt

## Report:
echo -e "\n#### haplotype_caller.sh: Done with script."
date
