#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Script adapted and modified from Poelstra et al. 2021, Systematic Biology (https://doi.org/10.1093/sysbio/syaa053)

## Software:
amas=/home/nibtve93/software/AMAS/amas/AMAS.py # (https://github.com/marekborowiec/AMAS)

## Command-line args:
fasta_dir=$1
locusstats=$2

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### 03e_locusstats.sh: Starting script."
echo -e "#### 03e_locusstats.sh: Input directory with FASTA files: $fasta_dir"
echo -e "#### 03e_locusstats.sh: Output file with statistics: $LOCUSTATS \n\n"

################################################################################
#### ESTIMATE STATISTICS FOR ALL LOCI  ####
################################################################################
echo -e "\n#### 03e_locusstats.sh: Estimate statistics for each locus ..."
for fasta in $fasta_dir/*
do
    echo -e "\n## 03e_locusstats.sh: Processing $fasta ..."
    fasta_id=$(basename $fasta)
    file_stats=$(dirname $locusstats)/tmp.$fasta_id.stats.txt

    python3 $amas summary -f fasta -d dna -i $fasta -o $file_stats

    grep -v "Alignment_name" $file_stats >> $locusstats
done

## Add header
echo -e "\n#### 03e_locusstats.sh: Adding header line ..."
header=$(head -n 1 $file_stats)
(echo $header && cat $locusstats) > $(dirname $locusstats)/tmp.txt && mv $(dirname $locusstats)/tmp.txt $locusstats

## Remove temporary files
echo -e "\n#### 03e_locusstats.sh: Removing temporary files ..."
find $(dirname $locusstats) -name 'tmp*txt' -delete

## Report:
echo -e "\n#### 03e_locusstats.sh: Done with script."
date


