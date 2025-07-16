#!/bin/sh
#SBATCH -p medium40

################################################################################
#### SET-UP ####
################################################################################
cd
source ~/.bashrc

## Software:
# iqtree2 needs to be included in $PATH (v2.2.0; https://iqtree.github.io/)

## Command-line args:
in_file=$1
prefix=$2
ufboot=$3

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### iqtree_bpp.sh: Starting script."
echo -e "#### iqtree_bpp.sh: Input file: $in_file"
echo -e "#### iqtree_bpp.sh: Prefix: $prefix"
echo -e "#### iqtree_bpp.sh: Number of ultrafast bootstrap replicates: $ufboot \n\n"

################################################################################
#### MAXIMUM LIKELIHOOD INFERENCE IN IQ-TREE 2 ####
################################################################################
iqtree2 -T AUTO -s $in_file --seqtype DNA -m GTR+G -nstop 200 -B $ufboot -wbt -bnni -alrt 1000 --prefix $prefix
