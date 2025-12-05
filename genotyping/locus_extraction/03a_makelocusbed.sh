#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################

## Command-line args:
scripts_dir=$1
set_id=$2
ind_file=$3
bed_dir=$4
locusbed_intermed=$5
min_elem_ovl=$6
min_elem_ovl_trim=$7
min_locus_size=$8
max_dist_within_ind=$9
max_dist_between_ind=$10
min_elem_size=$11
last_row=$12

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### 03a_makelocusbed.sh: Starting script."
echo -e "#### 03a_makelocusbed.sh: Directory with scripts: $scripts_dir"
echo -e "#### 03a_makelocusbed.sh: Set ID: $set_id"
echo -e "#### 03a_makelocusbed.sh: File with individuals: $ind_file"
echo -e "#### 03a_makelocusbed.sh: Directory with BED files: $bed_dir"
echo -e "#### 03a_makelocusbed.sh: Output BED file: $locusbed_intermed"
echo -e "#### 03a_makelocusbed.sh: Minimum element overlap for locus creation: $min_elem_ovl"
echo -e "#### 03a_makelocusbed.sh: Minimum element overlap for locus trimming: $min_elem_ovl_trim"
echo -e "#### 03a_makelocusbed.sh: Minimum locus size: $min_locus_size"
echo -e "#### 03a_makelocusbed.sh: Maximum distance within individuals: $max_dist_within_ind"
echo -e "#### 03a_makelocusbed.sh: Maximum distance between individuals: $max_dist_between_ind"
echo -e "#### 03a_makelocusbed.sh: Minimum locus size: $min_elem_size"
echo -e "#### 03a_makelocusbed.sh: Number of loci to process (all if 0): $last_row \n\n"

################################################################################
#### MAKE BED FILE WITH LOCUS COORDINATES ####
################################################################################
## Run GATK CallableLoci to produce BED file for sites that are (non-)callable for a single sample
echo -e "#### 03a_makelocusbed.sh: Running script to create BED file with locus coordinates ..."
Rscript $scripts_dir/03a_makelocusbed.R $set_id $ind_file $bed_dir $locusbed_intermed \
	$min_elem_ovl $min_elem_ovl_trim $min_locus_size $max_dist_within_ind $max_dist_between_ind $min_elem_size $last_row

## Report:
echo -e "\n#### 03a_makelocusbed.sh: Done with script."
date