#!/bin/bash
#SBATCH -p medium40
#SBATCH -t 48:00:00
#SBATCH --signal=B:12@1800

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# dmtcp needs to be included in $PATH (https://dmtcp.sourceforge.io/)
# bpp needs to be included in $PATH (v4.4.1; https://github.com/bpp/bpp)

## Command-line args:
check_in=$1
check_out=$2

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### bpp_checkpoint.sh: Starting script."
echo -e "#### bpp_checkpoint.sh: Previous checkpointing directory: $check_in"
echo -e "#### bpp_checkpoint.sh: Checkpointing directory: $check_out \n\n"

################################################################################
#### COALESCENT MODELLING IN BPP ####
################################################################################
trap 'echo -e "#### bpp_checkpoint.sh: Checkpointing ..."; date; dmtcp_command --bcheckpoint; echo -e "#### bpp_checkpoint.sh: Checkpointing done."; date; exit 12' 12

echo -e "#### bpp_checkpoint.sh: Coalescent modelling in BPP ..."
dmtcp_restart --ckptdir $check_out $check_in/ckpt_*.dmtcp &
wait