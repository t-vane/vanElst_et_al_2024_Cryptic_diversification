#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# BWA needs to be included in $PATH (v7.0.17; https://github.com/lh3/bwa)

## Command-line args:
mode=$1
nt=$2
index=$3
in_dir=$4
bam_dir=$5
in_file=$6
indv=$(sed -n "$SLURM_ARRAY_TASK_ID"p $in_file)
readgroup="@RG\tID:group1\tSM:$indv\tPL:illumina\tLB:lib1"

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### reference_mapping.sh: Starting script."
echo -e "#### reference_mapping.sh: Paired-end or single-end: $mode"
echo -e "#### reference_mapping.sh: Number of threads: $nt"
echo -e "#### reference_mapping.sh: Reference genome index: $index"
echo -e "#### reference_mapping.sh: Directory with trimmed reads: $in_dir"
echo -e "#### reference_mapping.sh: Directory for BAM files: $bam_dir"
echo -e "#### reference_mapping.sh: File with individuals: $in_file"
echo -e "#### reference_mapping.sh: Individual: $indv \n\n"

################################################################################
#### MAP TO REFERENCE GENOME ####
################################################################################
if [[ $mode == "PE" ]]
then
	echo -e "#### reference_mapping.sh: Reference mapping for paired-end individual $indv ...\n"
	bwa mem -aM -R $readgroup -t $nt $index $in_dir/$indv.trimmed.1.fq.gz $in_dir/$indv.trimmed.2.fq.gz | samtools view -b -h > $bam_dir/$indv.bam
elif [[ $mode == "SE" ]]
then
	echo -e "#### reference_mapping.sh: Reference mapping for single-end individual $indv ...\n"
	bwa mem -aM -R $readgroup -t $nt $index $in_dir/$indv.trimmed.1.fq.gz | samtools view -b -h > $bam_dir/$indv.bam
else
	echo -e "#### reference_mapping.sh: Invalid sequencing mode provided - only PE and SE allowed. ...\n" && exit 1
fi

## Report:
echo -e "\n#### reference_mapping.sh: Done with script."
date





