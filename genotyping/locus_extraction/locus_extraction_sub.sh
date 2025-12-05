################################################################################
#### LOCUS EXTRACTION WITH CUSTOM PIPELINE ####
################################################################################
## Pipeline adapted and modified from Poelstra et al. 2021, Systematic Biology (https://doi.org/10.1093/sysbio/syaa053)

scripts_dir=/home/nibtve93/scripts/locusExtraction

prefix=MicrocebusPhylogenomics
set_id=reducedSet
reference=$PWORK/mmur3/GCF_000165445.2_Mmur_3.0_genomic_reduced.fna
bam_dir=$PWORK/$prefix/bamLinks
vcf_dir=$PWORK/$prefix/gatk

# These are the VCF files with all individuals
vcf_altreff_full=$vcf_dir/allScaffolds.annot.SNP.vcf # Unfiltered SNP file
vcf_filt_mask_full=$vcf_dir/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.vcf # Filtered SNP file
# These will be the expected VCF files with selected individuals
vcf_altref=$vcf_dir/allScaffolds.annot.SNP.$set_id.vcf 
vcf_filt_mask=$vcf_dir/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.$set_id.vcf

out_dir=$PWORK/$prefix/gatk/locusExtraction/$set_id

fasta_dir=$out_dir/fasta
indfasta_dir=$fasta_dir/$set_id.byind
bed_dir=$out_dir/bed
bed_removed_sites=$bed_dir/$set_id.sitesinvcfremovedbyfilters.bed

locusfasta_dir_intermed=$fasta_dir/"$set_id"_bylocus_intermed/
locusfasta_dir_final=$fasta_dir/"$set_id"_bylocus_final

locusbed_intermed=$bed_dir/"$set_id"_loci_intermed.bed
locusbed_final=$bed_dir/"$set_id"_loci_all.bed
locuslist=$bed_dir/locuslist.txt
fasta_merged=$indfasta_dir/"$set_id"_merged.fasta

locusstats_intermed=$bed_dir/$set_id/"$set_id"_locusstats_all.txt
locusstats_final=$bed_dir/$set_id/"$set_id"_locusstats_filt.txt

mkdir -p $fasta_dir
mkdir -p $indfasta_dir
mkdir -p $bed_dir

mkdir $locusfasta_dir_intermed
mkdir $locusfasta_dir_final
mkdir $bed_dir/$set_id

mkdir -p $out_dir/logFiles/

ind_file=$out_dir/$set_id.txt # Contains individuals selected for divergence time estimation
no_inds=$(cat $ind_file | wc -l)

# Subset VCF files to selected individuals
bcftools view -S $ind_file -O v -o $vcf_altref $vcf_altref_full
bcftools view -S $ind_file $vcf_filt_mask_full > $vcf_filt_mask

#################################################################
#### 1 CREATE MASKED REFERENCE GENOME PER INDIVIDUAL ####
#################################################################
## Create BED file with masked sites from VCF
sbatch --account=nib00015 --output=$out_dir/logFiles/01_maskbed.$set_id.oe $scripts_dir/01_maskbed.sh $vcf_altref $vcf_filt_mask $bed_removed_sites $bed_dir

## Produce masked FASTA file per individual
suffix=auto
min_dp=2
sbatch --account=nib00015 --array=1-$no_inds --output=$out_dir/logFiles/02_process-inds.%A_%a.$set_id.oe $scripts_dir/02_process-inds.sh \
	$ind_file $vcf_altref $reference $bam_dir $suffix $min_dp $indfasta_dir $bed_dir $bed_removed_sites

#################################################################
#### 2 EXTRACT AND FILTER LOCI ACROSS INDIVIDUALS ####
#################################################################
## Make BED file with desired locus coordinates
min_elem_ovl=0.9 # Minimum element overlap for locus creation
min_elem_ovl_trim=0.8 # Minimum element overlap for locus trimming
min_locus_size=100 # Minimum locus size
max_dist_within_ind=10 # Maximum distance within individuals
max_dist_between_ind=0 # Maximum distance between individuals
min_elem_size=25 # Minimum locus size
last_row=0 # Number of loci to process (all if 0)
sbatch --account=nib00015 --output=$out_dir/logFiles/03a_makelocusbed.$set_id.oe $scripts_dir/03a_makelocusbed.sh $scripts_dir $set_id $ind_file $bed_dir $locusbed_intermed \
	$min_elem_ovl $min_elem_ovl_trim $min_locus_size $max_dist_within_ind $max_dist_between_ind $min_elem_size $last_row

## Intersection of BED file with loci with too high depth not necessary since these were already filtered out earlier (see genotype call filtering)
locusbed_final=$locusbed_intermed

## Get merged FASTA file with all individuals and loci
sbatch --account=nib00015 --output=$out_dir/logFiles/03c_mergedfasta.$set_id.oe $scripts_dir/03c_mergedfasta.sh $ind_file $indfasta_dir $locusbed_final $locuslist $fasta_merged

## Create by-locus FASTA files
nloci=$(cat $locuslist | wc -l)
sbatch --account=nib00015 --array=1-$nloci --output=$out_dir/logFiles/03d_locusfasta.%A_%a.$set_id.oe $scripts_dir/03d_locusfasta.sh $locuslist $locusfasta_dir_intermed $fasta_merged

## Estimate statistics for intermediate loci
sbatch --account=nib00015 --output=$out_dir/logFiles/03e_locusstats_intermed.$set_id.oe $scripts_dir/03e_locusstats.sh $locusfasta_dir_intermed $locusstats_intermed

## Filter loci for maximum proportion of missing data and minimum distance between loci
maxmiss=15 # Maximum percentage of missing data in percent
mindist=10000 # Minimum distance (bp) between loci
sbatch --account=nib00015 --output=$out_dir/logFiles/03f_filterloci.$set_id.oe $scripts_dir/03f_filterloci.sh $locusstats_intermed $locusfasta_dir_intermed $locusfasta_dir_final $maxmiss $mindist

## Estimate statistics for final loci
sbatch --account=nib00015 --output=$out_dir/logFiles/03e_locusstats_final.$set_id.oe $scripts_dir/03e_locusstats.sh $locusfasta_dir_final $locusstats_final

## Archive and remove intermediate locus files
for i in 03366 03367 03368 03369
do
tar -vcf $locusfasta_dir_intermed/${set_id}_loci_intermed_$i.tar $locusfasta_dir_intermed/NC_$i*fa && rm $locusfasta_dir_intermed/NC_$i*fa
done



