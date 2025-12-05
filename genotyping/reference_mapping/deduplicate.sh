#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# SAMtools needs to be included in $PATH (v1.11; http://www.htslib.org/)

## Command-line args:
bam_dir=$1
in_file=$2
minmapq=$3
indv=$(sed -n "$SLURM_ARRAY_TASK_ID"p $in_file)

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### deduplicate.sh: Starting script."
echo -e "#### deduplicate.sh: Directory for BAM files: $bam_dir"
echo -e "#### deduplicate.sh: File with individuals: $in_file"
echo -e "#### deduplicate.sh: Minimum mapping quality for filtering: $minmapq"
echo -e "#### deduplicate.sh: Individual: $indv \n\n"

################################################################################
#### DEDUPLICATE ####
################################################################################
echo -e "#### deduplicate.sh: Deduplication for paired-end individual $indv ...\n"
samtools collate --output-fmt BAM $bam_dir/$indv.MQ$minmapq.pp.bam -O | samtools fixmate - - -r -m -O BAM| samtools sort -m 15G -O BAM | samtools markdup - $bam_dir/$indv.MQ$minmapq.pp.dedup.bam -r -s -O BAM

## Report:
echo -e "\n#### deduplicate.sh: Done with script."
date
