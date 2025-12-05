#!/bin/bash

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# paup needs to be included in $PATH (v4.0a; https://paup.phylosolutions.com/)

## Command-line args:
in_file=$1
log_file=$2

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### svdq.sh: Starting script."
echo -e "#### svdq.sh: PAUP input file: $in_file"
echo -e "#### svdq.sh: Log file: $log_file \n\n"

################################################################################
#### PHYLOGENETIC INFERENCE WITH SVDquartets####
################################################################################
cd $(dirname $in_file)

echo -e "#### svdq.sh: Phylogenetic inference with SVDquartets... \n"
paup4a168_ubuntu64 -n $in_file $log_file

## Report:
echo -e "\n#### svdq.sh: Done with script."
date
