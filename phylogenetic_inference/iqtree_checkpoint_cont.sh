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
check_in=$1
check_out=$2

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### iqtree_checkpoint_cont.sh: Starting script."
echo -e "#### iqtree_checkpoint_cont.sh: Previous checkpointing directory: $check_in"
echo -e "#### iqtree_checkpoint_cont.sh: Checkpointing directory: $check_out \n\n"

################################################################################
#### Continue phylogenetic inference with IQ-TREE ####
################################################################################
trap 'echo -e "#### iqtree_checkpoint_cont.sh: Checkpointing ..."; date; dmtcp_command --bcheckpoint; echo -e "#### iqtree_checkpoint_cont.sh: Checkpointing done."; date; exit 12' 12

echo -e "#### iqtree_checkpoint_cont.sh: Coalescent modelling in G-PhoCS ..."
dmtcp_restart --ckptdir $check_out $check_in/ckpt_*.dmtcp &
wait

