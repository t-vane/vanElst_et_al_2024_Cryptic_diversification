#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# BWA needs to be included in $PATH (v7.0.17; https://github.com/lh3/bwa)
# SAMtools needs to be included in $PATH (v1.11; http://www.htslib.org/)
# gatk needs to be included in $PATH (v4.1.9.0; https://gatk.broadinstitute.org/hc/en-us)

## Command-line args:
reference=$1
bwa=$2
bwa_index=$3
samtools=$4
gatk=$5

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### indexing.sh: Starting script."
echo -e "#### indexing.sh: Reference genome: $reference"
echo -e "#### indexing.sh: Index with BWA: $bwa"
[[ $bwa ]] && echo -e "#### indexing.sh: Name for BWA index file: $bwa_index"
echo -e "#### indexing.sh: Index with SAMtools: $samtools"
echo -e "#### indexing.sh: Index with GATK (Picard): $gatk \n\n"

################################################################################
#### CREATE INDICES ####
################################################################################
cd $(dirname $reference)

## Index with BWA
[[ $bwa ]] && echo -e "#### indexing.sh: Creating index of $reference with BWA ...\n" && bwa index -p $bwa_index -a bwtsw $reference

## Index with SAMtools
[[ $samtools ]] && echo -e "#### indexing.sh: Creating index of $reference with SAMtools ...\n" && samtools faidx $reference

## Index with GATK (Picard)
if [[ $gatk ]]
then
	echo -e "#### indexing.sh: Creating index of $reference with GATK (Picard) ...\n"
	gatk CreateSequenceDictionary -R $reference -O $reference.dict
	oldname=${reference}.dict
	newname="${oldname/fna.dict/dict}"
	cp $oldname $newname
fi

## Report:
echo -e "\n#### indexing.sh: Done with script."
date
