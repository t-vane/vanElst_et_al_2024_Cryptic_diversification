################################################################################
#### SPECIES DELIMITATION THROUGH ANCESTRY INFERENCE ####
################################################################################
scripts_dir=/home/nibtve93/scripts/populationStructure

prefix=MicrocebusPhylogenomics

#################################################################
#### ANCESTRY INFERENCE WITH NGSADMIX ####
#################################################################

## Assign sample sets
sample_sets="Marn_Msp1 Mber_Mmyo_Mruf Mbon_Mdan_Mrav Mbor_Msim Mgan_Mman_Mmur Mger_Mjol_Mmaro Mjon_Mmac Mleh_Mmit Mmam_Mmar_Msam"

## Run NgsAdmix pipeline for each sample set
for set_id in $sample_sets
do
echo "... processing set $sample_set"

angsd_dir=$PWORK/$prefix/angsd/$set_id
beagle=$angsd_dir/$set_id.beagle.gz # Genotype likelihood file created in genotype_likelihoods_sub.sh
out_dir=$PWORK/$prefix/speciesDelimitation/genomics/admixture/set_id

mkdir -p $out_dir/logFiles

clusters=5 # Maximum number of clusters to assume in admixture analysis
repeats=10 # Number of independent runs
percentage="75/100" # Minimum percentage of represented individuals
minind=$(( ($(zcat $beagle | head -1 | wc -w)/3-1) * $percentage )) # Minimum number of represented individuals
nt=80

## Submit array job to infer individual ancestries for each number of clusters (ranging from 2 to $clusters), using $repeats repetitions 
for k in $(seq 2 $clusters)
do
	sbatch --array=1-$repeats --output=$out_dir/logFiles/ngsadmix.K$k.$set_id.%A_%a.oe $scripts_dir/ngsadmix.sh $nt $k $beagle $out_dir $minind $set_id
done

## Print likelihood values file
like_file=$out_dir/likevalues.$set_id.txt # File for likelihoods summary
rm $like_file; touch $like_file
for k in $(seq 2 $clusters); 
do
	for seed in $(seq 1 $repeats)
	do
		until [ -f $out_dir/$set_id.K$k.seed$seed.qopt ]
		do
			sleep 5m
		done
		
		sbatch --account=nib00015 --output=$out_dir/logFiles/print_likes.$set_id.oe $scripts_dir/best_likes.sh $out_dir/$set_id.K$k.seed$seed.log $like_file $k $seed
	done
done

done

