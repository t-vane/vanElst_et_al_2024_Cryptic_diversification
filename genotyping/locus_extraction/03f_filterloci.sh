#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################

## Command-line args:#
scripts_dir=$1
locusstats=$2
in_dir=$3
out_dir=$4
maxmiss=$5
mindist=$6

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### 03f_filterloci.sh: Starting script."
echo -e "#### 03f_filterloci.sh: Directory with scripts: $scripts_dir"
echo -e "#### 03f_filterloci.sh: File with locus statistics: $locusstats"
echo -e "#### 03f_filterloci.sh: Input FASTA directory: $in_dir"
echo -e "#### 03f_filterloci.sh: Output FASTA directory: $out_dir"
echo -e "#### 03f_filterloci.sh: Maximum percentage of missing data: $maxmiss"
echo -e "#### 03f_filterloci.sh: Minimum distance (bp) between loci: $mindist \n\n"

################################################################################
#### FILTER LOCI FOR MISSING DATA AND DISTANCE ####
################################################################################
echo -e "#### 03f_filterloci.sh: Running script to filter for missing data and distance ..."
Rscript $scripts_dir/03f_filterloci.R $locusstats $in_dir $out_dir $maxmiss$ $mindist

## Report:
echo -e "\n#### 03f_filterloci.sh: Done with script."
date