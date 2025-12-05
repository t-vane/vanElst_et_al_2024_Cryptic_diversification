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
gatk3=/home/nibtve93/software/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar (v3.8.1.0; https://gatk.broadinstitute.org/hc/en-us)

## Command-line args:
ind_file=$1
vcf_altref=$2
reference=$3
bam_dir=$4
suffix=$5
min_dp=$6
indfasta_dir=$7
bed_dir=$8
bed_removed_sites=$9
indv=$(sed -n "$SLURM_ARRAY_TASK_ID"p $ind_file)

bam=$bam_dir/$indv.*$suffix.bam

callable_summary=$bed_dir/$indv.callableloci.sumtable.txt
bed_out=$bed_dir/$indv.callablelocioutput.bed
bed_noncallable=$bed_dir/$indv.noncallable.bed
bed_callable=$bed_dir/$indv.callable.bed

fasta_altref=$indfasta_dir/$indv.altref.fasta
fasta_masked=$indfasta_dir/$indv.altrefmasked.fasta

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### 02_process-inds.sh: Starting script."
echo -e "#### 02_process-inds.sh: File with individuals: $ind_file"
echo -e "#### 02_process-inds.sh: Raw VCF: $vcf_altref"
echo -e "#### 02_process-inds.sh: Reference genome: $reference"
echo -e "#### 02_process-inds.sh: Directory with BAM files: $bam_dir"
echo -e "#### 02_process-inds.sh: BAM file suffix: $suffix"
echo -e "#### 02_process-inds.sh: Minimum depth for CallableLoci: $min_dp"
echo -e "#### 02_process-inds.sh: Directory for FASTA files: $indfasta_dir"
echo -e "#### 02_process-inds.sh: Directory with BED file: $bed_dir"
echo -e "#### 02_process-inds.sh: BED file with removed sites: $bed_removed_sites"
echo -e "#### 02_process-inds.sh: Individual: $indv"
echo -e "#### 02_process-inds.sh: BAM file: $bam"
echo -e "#### 02_process-inds.sh: CallableLoci summary: $callable_summary"
echo -e "#### 02_process-inds.sh: CallableLoci output BED all: $bed_out"
echo -e "#### 02_process-inds.sh: CallableLoci output BED non-callable: $bed_noncallable"
echo -e "#### 02_process-inds.sh: CallableLoci output BED callable: $bed_callable"
echo -e "#### 02_process-inds.sh: Raw FASTA genome: $fasta_altref"
echo -e "#### 02_process-inds.sh: Masked FASTA genome: $fasta_masked \n\n"

################################################################################
#### CREATE MASKED FASTA GENOME FOR INDIVIDUAL ####
################################################################################
## Index BAM file if not yet done
[[ ! -e $bam.bai ]] && echo -e "#### 02_process-inds.sh: Indexing bam ..." && samtools index $bam

## Run GATK CallableLoci to produce BED file for sites that are (non-)callable for a single sample
echo -e "#### 02_process-inds.sh: Running CallableLoci ..."
java -jar $gatk3 -T CallableLoci -R $reference -I $bam -summary $callable_summary --minDepth $min_dp -o $bed_out
# Separate into two files
grep -v "CALLABLE" $bed_out > $bed_noncallable
grep "CALLABLE" $bed_out > $bed_callable

## Produce whole-genome FASTA file for a single sample
echo -e "#### 02_process-inds.sh: Running FastaAlternateReferenceMaker ..."
java -jar $gatk3 -T FastaAlternateReferenceMaker -IUPAC $indv -R $reference -V $vcf_altref -o $fasta_altref
# Edit FASTA headers
sed -i -e 's/:1//g' -e 's/>[0-9]* />/g' $fasta_altref
## Count bases
n_acgt=$(grep -Eo "A|C|G|T" $fasta_altref | wc -l)
n_ambig=$(grep -Eo "M|R|W|S|Y|K" $fasta_altref | wc -l)
echo -e "#### 02_process-inds.sh: Number of A/C/G/T in fasta_altref: $n_acgt"
echo -e "#### 02_process-inds.sh: Number of heterozygous sites (counted by ambiguity code) in fasta_altref: $n_ambig"

## Mask sites identified as non-callable or removed during VCF filtering in whole-genome FASTA
echo -e "#### 02_process-inds.sh: Running BEDtools maskfasta ..."
# Mask non-callable sites
echo -e "## 02_process-inds.sh: Non-callable sites ..."
bedtools maskfasta -fi $fasta_altref -bed $bed_noncallable-fo $fasta_masked.intermed.fasta
# Mask removed sites
echo -e "## 02_process-inds.sh: Callable sites ..."
bedtools maskfasta -fi $fasta_masked.intermed.fasta -bed $bed_removed_sites -fo $fasta_masked
# Counting Ns in FASTA files:
ncount_fasta_altref=$(grep -Fo "N" $fasta_altref | wc -l)
ncount_fasta_masked_intermed=$(grep -Fo "N" $fasta_masked.intermed.fasta | wc -l)
ncount_fasta_masked=$(grep -Fo "N" $fasta_masked | wc -l)
echo -e "#### 02_process-inds.sh: Number of Ns in fasta_altref: $ncount_fasta_altref"
echo -e "#### 02_process-inds.sh: Number of Ns in FASTA_MASKED_INTERMED: $ncount_fasta_masked_intermed"
echo -e "#### 02_process-inds.sh: Number of Ns in fasta_masked: $ncount_fasta_masked"

## Report:
echo -e "\n#### 02_process-inds.sh: Done with script."
date
