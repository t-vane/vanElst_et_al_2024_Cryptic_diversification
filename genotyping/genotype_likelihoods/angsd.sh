#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# angsd needs to be included in $PATH (v0.934; http://www.popgen.dk/angsd/index.php/ANGSD)

## Command-line args:
nt=$1
reference=$2
bamlist=$3
todo=$4
filters=$5
out=$6

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### angsd.sh: Starting script."
echo -e "#### angsd.sh: Number of threads: $nt"
echo -e "#### angsd.sh: Reference genome: $reference"
echo -e "#### angsd.sh: List of BAM files: $bamlist"
echo -e "#### angsd.sh: Inferences to do: $todo"
echo -e "#### angsd.sh: Filters to apply: $filters"
echo -e "#### angsd.sh: Output prefix: $out \n\n"

################################################################################
#### RUN ANGSD WITH SPECIFIED FILTERS AND INFERENCES ####
################################################################################
echo -e "#### angsd.sh: Running angsd ...\n"
angsd -nThreads $NC -ref $reference -bam $bamlist $todo $filters -out $out 

echo -e "\n#### angsd.sh: Done with script."
date

