#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Script adapted and modified from Poelstra et al. 2021, Systematic Biology (https://doi.org/10.1093/sysbio/syaa053)

## Software:
# BEDtools needs to be included in $PATH (v2.30.0; https://bedtools.readthedocs.io/en/latest/)
# SAMtools needs to be included in $PATH (v1.11; http://www.htslib.org/)

## Command-line args:
ind_file=$1
indfasta_dir=$2
locusbed_final=$$3
locuslist=$4
fasta_merged=$5

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### 03c_mergedfasta.sh: Starting script."
echo -e "#### 03c_mergedfasta.sh: File with individuals: $ind_file"
echo -e "#### 03c_mergedfasta.sh: Directory for FASTA files: $indfasta_dir"
echo -e "#### 03c_mergedfasta.sh: Final locus BED file: $locusbed_final"
echo -e "#### 03c_mergedfasta.sh: List of loci to be created: $locuslist"
echo -e "#### 03c_mergedfasta.sh: Merged output FASTA file: $fasta_merged \n\n"

################################################################################
#### CREATE MERGED FASTA FILE WITH ALL INDIVIDUALS AND LOCI ####
################################################################################
## Extract loci in locus BED file from masked FASTA for each individual
echo -e "#### 03c_mergedfasta.sh: Extracting loci from locus BED file from masked FASTA for each individual ..."
while read -r indv
do
	echo -e "## 03c_mergedfasta.sh: Processing $indv ..."
    fasta_in=$indfasta_dir/$indv.altrefmasked.fasta
    fasta_out=$indfasta_dir/${indv}_allloci.fasta
    bedtools getfasta -fi $fasta_in -bed $locusbed_final > $fasta_out
done < $ind_file

## Create list of loci
echo -e "#### 03c_mergedfasta.sh: Creating list of loci ..."
FASTA_1=$(find $indfasta_dir/*allloci.fasta | head -1)
grep ">" $FASTA_1 | sed 's/>//' > $locuslist

## Merge by-individual FASTA files
echo -e "#### 03c_mergedfasta.sh: Merging loci of individuals ..."
while read -r indv
do
	echo -e "## 03c_mergedfasta.sh: Processing $indv ..."
    fasta=$indfasta_dir/${indv}_allloci.fasta
    # Replace ":" by "," for compatibility with SAMtools faidx:
    sed "s/>/>${ind}__/g" $fasta | sed 's/:/,/g' >> $fasta_merged
done < $ind_file

## Index merged FASTA file with SAMtools faidx
echo -e "#### 03c_mergedfasta.sh: Indexing merged FASTA file with SAMtools faidx ..."
samtools faidx $fasta_merged

## Report:
echo -e "\n#### 03c_mergedfasta.sh: Done with script."
date

