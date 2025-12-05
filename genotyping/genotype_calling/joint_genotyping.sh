#!/bin/bash
#SBATCH -p medium40

set -euo pipefail

################################################################################
#### SET-UP ####
################################################################################
## Software:
# gatk needs to be included in $PATH (v4.1.9.0; https://gatk.broadinstitute.org/hc/en-us)

## Command-line args:
nt=$1
reference=$2
ind_file=$3
region_file=$4
gvcf_dir=$5
db_dir=$6
vcf_scaffold_dir=$7
tmp_dir=$8

scaffold_name=$(sed -n "$SLURM_ARRAY_TASK_ID"p $region_file | cut -f 1)
scaffold_start=$(sed -n "$SLURM_ARRAY_TASK_ID"p $region_file | cut -f 2)
scaffold_end=$(sed -n "$SLURM_ARRAY_TASK_ID"p $region_file | cut -f 3)

inds_command=$(for ind in `cat $ind_file`; do printf " --variant $gvcf_dir/$ind.rawvariants.g.vcf.gz"; done)

## Activate conda environment
conda activate java

## Report:
echo -e "\n\n###################################################################"
date
echo -e "#### joint_genotyping.sh: Starting script."
echo -e "#### joint_genotyping.sh: Number of threads: $nt"
echo -e "#### joint_genotyping.sh: Reference genome: $reference"
echo -e "#### joint_genotyping.sh: List with individuals: $ind_file"
echo -e "#### joint_genotyping.sh: List with regions and coordinates: $region_file"
echo -e "#### joint_genotyping.sh: Directory with GVCF files: $gvcf_dir"
echo -e "#### joint_genotyping.sh: Database directory: $db_dir"
echo -e "#### joint_genotyping.sh: Output directory for per-scaffold VCF files: $vcf_scaffold_dir"
echo -e "#### joint_genotyping.sh: Temporary directory for running GenotypeGVCFs: $tmp_dir"
echo -e "#### joint_genotyping.sh: Scaffold name: $scaffold_name"
echo -e "#### joint_genotyping.sh: Scaffold start coordinate: $scaffold_start"
echo -e "#### joint_genotyping.sh: Scaffold end coordinate: $scaffold_end \n\n"

################################################################################
#### CONDUCT JOINT GENOTYPING ####
################################################################################
echo -e "#### joint_genotyping.sh: Creating GenomicsDB data storage for scaffold $scaffold_name ...\n"
gatk GenomicsDBImport $inds_command --genomicsdb-shared-posixfs-optimizations --genomicsdb-workspace-path $db_dir/$scaffold_name --batch-size 0 -L $scaffold_name:$scaffold_start-$scaffold_end --reader-threads $nt --interval-padding 100

echo -e "#### joint_genotyping.sh: Conducting joint genotyping for scaffold $scaffold_name ...\n"
gatk GenotypeGVCFs -R $reference -V "gendb://$db_dir/$scaffold_name" -O $vcf_scaffold_dir/$scaffold_name.vcf.gz --tmp-dir $tmp_dir --genomicsdb-shared-posixfs-optimizations

## Report:
echo -e "\n#### joint_genotyping.sh: Done with script."
date

