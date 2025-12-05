#!/bin/bash
#SBATCH -p medium40
#SBATCH -t 48:00:00
#SBATCH --signal=B:12@1800

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
#dmtcp needs to be included in $PATH (https://dmtcp.sourceforge.io/)
#IQ-TREE needs to be included in $PATH (http://www.iqtree.org/)

module load gcc/8.2.0
export CC=/sw/compiler/gcc/8.2.0/skl/bin/gcc
export CXX=/sw/compiler/gcc/8.2.0/skl/bin/g++

## Command-line args:
in_file=$1
prefix=$2
ufboot=$3
check_out=$4

mkdir -p $check_out

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### gphocs_checkpoint.sh: Starting script."
echo -e "#### gphocs_checkpoint.sh: Input alignment: $in_file"
echo -e "#### gphocs_checkpoint.sh: Output prefix: $prefix"
echo -e "#### gphocs_checkpoint.sh: Number of ultrafast bootstrap replicates: $ufboot"
echo -e "#### gphocs_checkpoint.sh: Checkpointing directory: $check_out \n\n"

################################################################################
#### Phylogenetic inference with IQ-TREE ####
################################################################################
trap 'echo -e "#### iqtree_checkpoint.sh: Checkpointing ..."; date; dmtcp_command --bcheckpoint; echo -e "#### iqtree_checkpoint.sh: Checkpointing done."; date; exit 12' 12

echo -e "#### iqtree_checkpoint.sh: Phylogenetic inference with IQ-TREE ..."
dmtcp_launch --ckptdir $check_out iqtree2dyn -T AUTO -s $in_file --seqtype DNA -m GTR+G+ASC -nstop 200 -B $ufboot -wbt -bnni -alrt 1000 --prefix $prefix &
wait



