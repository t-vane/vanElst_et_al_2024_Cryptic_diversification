################################################################################
#### GENOTYPE CALLING WITH GATK####
################################################################################
scripts_dir=/home/nibtve93/scripts/gatk

prefix=MicrocebusPhylogenomics
vcf_dir=$PWORK/$prefix/gatk/
bam_dir=$PWORK/bamFiles/$prefix
ind_file=$vcf_dir/individuals.txt # Contains individual IDs in list format (".bam" will be appended in haplotypeCaller.sh)
reference_dir=$PWORK/mmur3
reference=$reference_dir/GCF_000165445.2_Mmur_3.0_genomic.fna # Reference genome in fasta format
region_file=$reference_dir/regionFileAutosomes_modified.bed # Regions (i.e., scaffolds) for which genotyping should be conducted jointly; should start at 1 and not at 0

gvcf_dir=$vcf_dir/gvcfFiles
db_dir=$vcf_dir/DBs
vcf_scaffold_dir=$vcf_dir/vcfFilesScaffolds
tmp_dir=$vcf_dir/tmp

mkdir -p $vcf_dir/logFiles
mkdir -p $tmp_dir
mkdir -p $gvcf_dir
mkdir -p $db_dir
mkdir -p $vcf_scaffold_dir

## Variant discovery with haplotype caller
nt=40
mem=100
no_inds=$(cat $ind_file | wc -l)
suffix=auto.bam
sbatch --array=1-$no_inds -c $nt --account=nib00015 --output=$vcf_dir/logFiles/haplotype_caller.%A_%a.oe $scripts_dir/haplotype_caller.sh $nt $mem $reference $ind_file $bam_dir $gvcf_dir $suffix

## Joint genotyping per scaffold
no_regions=$(cat $region_file | wc -l)
sbatch --array=1-$no_regions -c $nt --account=nib00015 --output=$vcf_dir/logFiles/joint_genotyping.%A_%a.oe $scripts_dir/joint_genotyping.sh $nt $reference $ind_file $region_file $gvcf_dir $db_dir $vcf_scaffold_dir $tmp_dir

## Merge per-scaffold VCFs
ls $vcf_scaffold_dir/*vcf.gz > $vcf_scaffold_dir/allScaffolds.vcflist
sbatch --account=nib00015 --output=$vcf_dir/logFiles/merge_vcfs.oe $scripts_dir/merge_vcfs.sh $vcf_scaffold_dir/allScaffolds.vcflist $vcf_dir/allScaffolds.vcf.gz
gzip -d $vcf_dir/allScaffolds.vcf.gz

################################################################################
#### VCF FILTERING ####
################################################################################
scripts_dir=/home/nibtve93/scripts/vcfFiltering

## Annotate with number of samples with data (NS INFO field) and remove indels and invariant sites
sbatch --account=nib00015 --output=$vcf_dir/logFiles/01_filter_indels_invariants.$prefix.oe $scripts_dir/01_filter_indels_invariants.sh $reference $vcf_dir/allScaffolds.vcf.gz $vcf_dir/allScaffolds.annot.SNP.vcf

## Apply filters based on sequencing depth/coverage, minor allele count and samples with data
in_file_all=$vcf_dir/inds_all.txt # List of all samples (without file extensions) that shall be included for global genotype likelihood estimation
no_inds_all=$(cat $in_file_all | wc -l)

mkdir -p $vcf_dir/bamHits

# Create bamHits file with locus coverages for each individual
for i in pe se
do
	in_file=$vcf_dir/inds_$i.txt # List of samples (without file extensions) that shall be included for global genotype likelihood estimation, separated by sequencing mode (pe or SE)
	no_inds=$(cat $in_file | wc -l)

	sbatch --array=1-$no_inds --output=$vcf_dir/logFiles/bamHits.$prefix.$i.%A_%a.oe $scripts_dir/coverage.sh $i $in_file $bam_dir $angsd_dir/bamHits
done

# Estimate and plot coverage distributions for each individual
sbatch --output=$vcf_dir/logFiles/cov_plot.$prefix.oe $scripts_dir/cov_plot.sh $scripts_dir $in_file_all $vcf_dir/bamHits $prefix.gc

# Set thresholds
mindepthind=$(cat $vcf_dir/bamHits/statistics/$prefix.gc.minmax.txt | cut -d " " -f2 | sort -n | head -1) # Minimum depth per individual
maxdepthind=$(cat $vcf_dir/bamHits/statistics/$prefix.gc.minmax.txt | cut -d " " -f3 | sort -n | tail -1) # Maximum depth per individual
gmin=10 # Minimum depth across all individuals
gmax=$(cat $vcf_dir/bamHits/statistics/$prefix.gc.minmax.txt | cut -d " " -f3 | paste -sd+ | bc) # Maximum depth across all individuals
minind=3 # Minimum number of represented individuals
mac=2 # Minor allele count

# Run filtering script
sbatch --account=nib00015 --output=$vcf_dir/logFiles/02_filter_dp.$prefix.oe $scripts_dir/02_filter_dp.sh $mindepthind $maxdepthind $gmin $gmax $minind $mac $vcf_dir/allScaffolds.annot.SNP.vcf $vcf_dir/allScaffolds.annot.SNP.minInd.DP.mac.vcf 

## The following part of the pipeline has been adapted and modified from Poelstra et al. 2021, Systematic Biology (https://doi.org/10.1093/sysbio/syaa053)

## Annotate with INFO fields FisherStrand, RMSMappingQuality, MappingQualityRankSumTest, ReadPosRankSumTest and AlleleBalance
suffix=auto
sbatch --account=nib00015 --output=$vcf_dir/logFiles/03_annot_gatk.$prefix.oe $scripts_dir/03_annot_gatk.sh \
	$vcf_dir/allScaffolds.annot.SNP.minInd.DP.mac.vcf  $vcf_dir/allScaffolds.annot.SNP.minInd.DP.mac.GATK.vcf  $bam_dir $reference $suffix

## Retain only bi-allelic sites and filter for INFO fields FisherStrand, RMSMappingQuality, MappingQualityRankSumTest, ReadPosRankSumTest and AlleleBalance
sbatch --wait --account=nib00015 --output=$vcf_dir/logFiles/04_filter_gatk.$prefix.oe $scripts_dir/04_filter_gatk.sh \
	$vcf_dir/allScaffolds.annot.SNP.minInd.DP.mac.GATK.vcf $vcf_dir/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-soft.vcf $vcf_dir/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.vcf $reference

## Final round of filtering to create VCF files with varying amounts of missing data across genotypes
filter_inds=FALSE # Boolean specifying whether to filter for missingness across individuals
maxmiss_ind=1 # Maxmimum missingness across individuals
for i in 0.05 0.25 0.5 0.75 0.95
do
maxmiss_geno=$(lua -e "print(1 - $i)") # Maximum missingness is inverted for genotypes, with a value of 1 indicating no missing data (this is how the '--max-missing' option in VCFtools works).
sbatch --account=nib00015 --output=$vcf_dir/logFiles/05_filter_missing.$prefix.oe $scripts_dir/05_filter_missing.sh \
	$vcf_dir/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.vcf $vcf_dir/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$i.vcf $maxmiss_geno $filter_inds $maxmiss_ind
done

