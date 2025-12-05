################################################################################
#### REFERENCE MAPPING AND FILTERING ####
################################################################################
scripts_dir=/home/nibtve93/scripts/referenceMapping

prefix=MicrocebusPhylogenomics
reference_dir=$PWORK/references/mmur3
reference=$reference_dir/GCF_000165445.2_Mmur_3.0_genomic.fna # Reference genome in fasta format

#################################################################
#### 0 INDEX REFERENCE GENOME IF NOT DONE YET ####
#################################################################
mkdir -p $reference_dir/logFiles

bwa=true # Boolean specifying whether to create index with BWA
bwa_index=mmur3 # Output name of bwa index
samtools=true # Boolean specifying whether to create index with SAMtools
gatk=true # Boolean specifying whether to create index with GATK (Picard)

## Submit indexing script
sbatch --output=$reference_dir/logFiles/indexing.$prefix.oe $scripts_dir/indexing.sh $reference $bwa $bwa_index $samtools $gatk

#################################################################
#### 1 ALIGN TRIMMED READS TO REFERENCE GENOME AND FILTER ####
#################################################################
nt=6

in_dir=$PWORK/trimmedReads/$prefix # The following scripts assume that trimmed read files are named *.trimmed.1.fq.gz and *.trimmed.2.fq.gz 
out_dir=$PWORK/bamFiles/$prefix

mkdir -p $out_dir/logFiles

## Submit scripts for reference mapping, sorting, filtering, extraction of genomic regions and header cleaning
for i in pe se
do
	in_file=$out_dir/map_$i.txt # List of samples (without file extensions) for which reference alignment shall be conducted
	no_inds=$(cat $in_file | wc -l)
	minmapq=20
	
	# Reference mapping, filtering and extraction of genomic regions
	sbatch --job-name=map_filter_pip --array=1-$no_inds --output=$out_dir/logFiles/reference_mapping_$i.%A_%a.$prefix.oe $scripts_dir/reference_mapping.sh $i $nt $reference_dir/$bwa_index $in_dir $out_dir $in_file
	
	# Sort and quality filter
	sbatch --job-name=map_filter_pip --dependency=singleton --array=1-$no_inds --output=$out_dir/logFiles/quality_filter_$i.%A_%a.$prefix.oe $scripts_dir/quality_filter.sh $i $nt $out_dir $in_file $minmapq
	
	# Deduplicate 
	[[ $i == pe ]] && sbatch --job-name=map_filter_pip --dependency=singleton --array=1-$no_inds --output=$out_dir/logFiles/deduplicate_$i.%A_%a.$prefix.oe $scripts_dir/deduplicate.sh $out_dir $in_file $minmapq
	
	# Extract genomic regions
	bed=$out_dir/regionFile_autosomes.bed # BED file with genomic regions that shall be extracted
	exclude="NW_|NC_028718.1|NC_033692.1" # String of chromosomes that are no longer represented (separator: "|")
	suffix=auto # Suffix for naming of final BAM files
	sbatch --job-name=map_filter_pip --dependency=singleton --array=1-$no_inds --output=$out_dir/logFiles/extract_regions_$i.%A_%a.$prefix.oe $scripts_dir/extract_regions.sh $i $out_dir $in_file $minmapq $bed "$exclude" $suffix
done
