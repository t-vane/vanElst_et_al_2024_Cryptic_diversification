#!/bin/sh
#SBATCH -p medium40

################################################################################
#### SET-UP ####
################################################################################
cd
source /home/nibtve93/.bashrc

## Software:
# iqtree2 needs to be included in $PATH (v2.2.0; https://iqtree.github.io/)

## Command-line args:
fasta_dir=$1
locus_dir=$2
prefix=$3

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### prepare_loci.sh: Starting script."
echo -e "#### prepare_loci.sh: Input FASTA directory: $in_file"
echo -e "#### prepare_loci.sh: Output locus directory: $locus_dir"
echo -e "#### prepare_loci.sh: Prefix: $prefix \n\n"

################################################################################
#### LOCUS PROCESSING ####
################################################################################
for fasta_id_long in $fasta_dir/*fa
do

fasta_id=$(basename $fasta_id_long)
basefile=$locus_dir/$prefix.$fasta_id

## Remove locus names from fasta headers and prefix sample name with "^"
sed 's/__.*//' $fasta_dir/$fasta_id | sed -e 's/.*_\(.*_A[01]\)_.*/>\1/' | sed -e 's/>/>^/' > $basefile.tmp1

## Transform to one-sample-per-line format (tmp 2); remove second and third columns
faidx --transform transposed $basefile.tmp1 | tee $basefile.tmp2 | cut -f 1,4 > $basefile.tmp3

## Estimate number of bases
seq_len=$(cut -f 3 $basefile.tmp2 | head -n 1)

## Estimate number of individuals
n_seqs=$(cat $basefile.tmp2 | wc -l)

## Print locus files
if [ $n_seqs != 0 ]
then
if [ $seq_len != 0 ]
then
echo "$n_seqs $seq_len" > $basefile.locus
cat $basefile.tmp3 >> $basefile.locus
else
echo "#### prepare_loci.sh: Skipping locus $fasta_id (sequence length = 0)..."
fi
else
echo "#### prepare_loci.sh: Skipping locus $fasta_id (no sequences)..."
fi
done

rm -f $locus_dir/*tmp*
rm -f $locus_dir/*fai