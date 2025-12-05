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
bam_dir=$2
in_file=$3
minmapq=$4
bed=$5
exclude=$6
suffix=$7
indv=$(sed -n "$SLURM_ARRAY_TASK_ID"p $in_file)

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### extract_regions.sh: Starting script."
echo -e "#### extract_regions.sh: Paired-end or single-end: $mode"
echo -e "#### extract_regions.sh: Directory for BAM files: $bam_dir"
echo -e "#### extract_regions.sh: File with individuals: $in_file"
echo -e "#### extract_regions.sh: Minimum mapping quality for filtering: $minmapq"
echo -e "#### extract_regions.sh: BED file with information on genomic regions: $bed"
echo -e "#### extract_regions.sh: String with chromosomes to exclude from header: $exclude"
echo -e "#### extract_regions.sh: Suffix for final BAM files: $suffix"
echo -e "#### extract_regions.sh: Individual: $indv \n\n"

################################################################################
#### EXTRACT SPECIFIC REGIONS FROM BED FILE ####
################################################################################
if [[ $mode == "PE" ]]
then
	echo -e "#### extract_regions.sh: Extraction of genomic regions provided in $bed for paired-end individual $indv ...\n"
	samtools view -b -L $bed $bam_dir/$indv.MQ$minmapq.pp.dedup.bam > $bam_dir/$indv.MQ$minmapq.pp.dedup.$suffix.bam
	
	echo -e "#### extract_regions.sh: Removing chromosomes $exclude from header for paired-end individual $indv ...\n"
	samtools view -H $bam_dir/$indv.MQ$minmapq.pp.dedup.$suffix.bam | egrep -v $exclude > $bam_dir/$indv.fixedhead
	cat $bam_dir/$indv.fixedhead <(samtools view $bam_dir/$indv.MQ$minmapq.pp.dedup.$suffix.bam) | samtools view -bo $bam_dir/$indv.MQ$minmapq.pp.dedup.$suffix.bam
	rm $bam_dir/$indv.fixedhead
elif [[ $mode == "SE" ]]
	echo -e "#### extract_regions.sh: Extraction of genomic regions provided in $bed for single-end individual $indv ...\n"
	samtools view -b -L $bed $bam_dir/$indv.MQ$minmapq.bam > $bam_dir/$indv.MQ$minmapq.$suffix.bam
	
	echo -e "#### extract_regions.sh: Removing chromosomes $exclude from header for single-end individual $indv ...\n"
	samtools view -H $bam_dir/$indv.MQ$minmapq.$suffix.bam | egrep -v $exclude > $bam_dir/$indv.fixedhead
	cat $bam_dir/$indv.fixedhead <(samtools view $bam_dir/$indv.MQ$minmapq.$suffix.bam) | samtools view -bo $bam_dir/$indv.MQ$minmapq.$suffix.bam
	rm $bam_dir/$indv.fixedhead
else
	echo -e "#### extract_regions.sh: Invalid sequencing mode provided - only PE and SE allowed. ...\n" && exit 1
fi

## Report:
echo -e "\n#### extract_regions.sh: Done with script."
date
