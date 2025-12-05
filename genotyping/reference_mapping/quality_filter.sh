#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# SAMtools needs to be included in $PATH (v1.11; http://www.htslib.org/)

## Command-line args:
mode=$1
nt=$2
bam_dir=$3
in_file=$4
minmapq=$5
indv=$(sed -n "$SLURM_ARRAY_TASK_ID"p $in_file)

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### quality_filter.sh: Starting script."
echo -e "#### quality_filter.sh: Paired-end or single-end: $mode"
echo -e "#### quality_filter.sh: Number of threads: $nt"
echo -e "#### quality_filter.sh: Directory for BAM files: $bam_dir"
echo -e "#### quality_filter.sh: File with individuals: $in_file"
echo -e "#### quality_filter.sh: Minimum mapping quality for filtering: $minmapq"
echo -e "#### quality_filter.sh: Individual: $indv \n\n"

################################################################################
#### SORT AND FILTER FOR MINIMUM MAPPING QUALITY AND PROPER PAIRING (IF PE) ####
################################################################################
if [[ $mode == "PE" ]]
then
	echo -e "#### quality_filter.sh: Minimum mapping quality and proper-pair filtering and sorting for paired-end individual $indv ...\n"
	samtools view -bhu -q $minmapq -f 0x2 -@ $nt $bam_dir/$indv.bam | samtools sort -@ $nt -m 15G -O bam > $bam_dir/$indv.MQ$minmapq.pp.bam
elif [[ $mode == "SE" ]]
	echo -e "#### quality_filter.sh: Minimum mapping quality filtering and sorting for single-end individual $indv ...\n"
	samtools view -bhu -q $minmapq -@ $nt $bam_dir/$indv.bam | samtools sort -@ $nt -m 15G -O bam > $bam_dir/$indv.MQ$minmapq.bam
else
	echo -e "#### quality_filter.sh: Invalid sequencing mode provided - only PE and SE allowed. ...\n" && exit 1
fi

## Report:
echo -e "\n#### quality_filter.sh: Done with script."
date