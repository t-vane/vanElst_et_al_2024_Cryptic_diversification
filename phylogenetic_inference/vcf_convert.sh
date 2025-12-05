#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# VCFtools needs to be included in $PATH (v0.1.17; https://vcftools.github.io/index.html)
amas=/home/nibtve93/software/AMAS/amas/AMAS.py # (https://github.com/marekborowiec/AMAS)
ascbias=/home/nibtve93/software/raxml_ascbias/ascbias.py # (https://github.com/btmartin721/raxml_ascbias/blob/master/ascbias.py)

## Command-line args:
scripts_dir=$1
vcf_in=$2
out_dir=$3
format=$4

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### vcf_convert.sh: Starting script."
echo -e "#### vcf_convert.sh: Directory with vcf_tab_to_fasta_alignment script: $scripts_dir"
echo -e "#### vcf_convert.sh: Input VCF file: $vcf_in"
echo -e "#### vcf_convert.sh: Output directory: $out_dir"
echo -e "#### vcf_convert.sh: Output format of alignment: $format \n\n"

################################################################################
#### CONVERT VCF TO PHYLIP FORMAT AND CALCULATE STATISTICS ####
################################################################################
PREFIX=$(basename $vcf_in .vcf)

if [[ $format == phylip ]]
then
	suffix=phy
elif [[ $format == nexus ]]
then
	suffix=nex
else
	echo -e "#### vcf_convert.sh: Invalid format provided - only phylip and nexus accepted. ...\n" && exit 1
fi

echo -e "#### vcf_convert.sh: Creating TAB file from VCF ...\n"
vcf-to-tab < $vcf_in > $out_dir/$prefix.tab

echo -e "#### vcf_convert.sh: Converting TAB file to FASTA format ...\n"
perl $scripts_dir/vcf_tab_to_fasta_alignment_TVE.pl -i $out_dir/$prefix.tab > $out_dir/$prefix.fasta

echo -e "#### vcf_convert.sh: Converting FASTA file to $format format ...\n"
python $amas convert -i $out_dir/$prefix.fasta -f fasta -u $format -d dna
mv $out_dir/$prefix.fasta-out.$suffix $out_dir/$prefix.$suffix

echo -e "#### vcf_convert.sh: Removing invariant sites ...\n"
python $ascbias -p $out_dir/$prefix.$suffix -o $out_dir/$prefix.noinv.$suffix

echo -e "#### vcf_convert.sh: Calculating alignment statistics ...\n"
python $amas summary -f $format -d dna -i $out_dir/$prefix.noinv.$suffix -o $out_dir/$prefix.noinv.$suffix.summary

echo -e "#### vcf_convert.sh: Writing partitions file ...\n"
align_length=$(sed -n 2p $out_dir/$prefix.noinv.$suffix.summary | cut -f3)
echo "[asc~$out_dir/$prefix.noinv.$suffix.stamatakis], ASC_DNA, p1=1-$align_length" > $out_dir/$prefix.noinv.$suffix.partitions

## Report:
echo -e "\n#### vcf_convert.sh: Done with script."
date


