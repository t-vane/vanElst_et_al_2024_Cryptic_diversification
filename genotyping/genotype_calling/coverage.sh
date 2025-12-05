#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# samtools needs to be included in $PATH (v1.11; http://www.htslib.org/)

## Command-line args:
mode=$1
in_file=$2
bam_dir=$3
out_dir=$4
indv=$(sed -n "$SLURM_ARRAY_TASK_ID"p $in_file)

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### coverage.sh: Starting script."
echo -e "#### coverage.sh: Sequencing mode: $mode"
echo -e "#### coverage.sh: File with individuals: $in_file"
echo -e "#### coverage.sh: Directory with BAM files: $bam_dir"
echo -e "#### coverage.sh: Output directory: $out_dir"
echo -e "#### coverage.sh: Individual: $indv \n\n"

################################################################################
#### ESTIMATE NUMBER OF SITES AND COVERAGE ####
################################################################################
if [[ $mode == "PE" ]]
then
	echo -e "#### coverage.sh: Estimating coverage for each SbfI forward strand site of $indv ...\n"
	samtools view -f 0x40 $bam_dir/$indv.auto.bam | awk '$10 ~ /^TGCAGG/' | cut -f 3,4 | tr '\t' '_' | sort | uniq -c | sort -k1 -nr | sed -E 's/^ *//; s/ /\t/' > $out_dir/$indv.sbf1.f1.bamhits
elif [[ $mode == "SE" ]]
then
	echo -e "#### coverage.sh: Estimating coverage for each SbfI forward strand site of $indv ...\n"
	samtools view $bam_dir/$indv.auto.bam | awk '$10 ~ /^TGCAGG/' | cut -f 3,4 | tr '\t' '_' | sort | uniq -c | sort -k1 -nr | sed -E 's/^ *//; s/ /\t/' > $out_dir/$indv.sbf1.f1.bamhits
else
	echo -e "#### coverage.sh: Invalid sequencing mode provided - only PE and SE allowed. ...\n" && exit 1
fi

## Report:
nsites=$(wc -l $out_dir/$indv.sbf1.f1.bamhits)
nreads=$(cat $out_dir/$indv.sbf1.f1.bamhits | cut -f 1 | paste -sd+ | bc)
mean_cov=$(( $nsites / $nreads ))

echo -e "#### coverage.sh: Number of SbfI sites recovered for $indv: $nsites"
echo -e "#### coverage.sh: Total number of reads across these sites: $nreads"
echo -e "#### coverage.sh: Mean coverage across these sites: $mean_cov"

echo -e "\n#### coverage.sh: Done with script."
date



