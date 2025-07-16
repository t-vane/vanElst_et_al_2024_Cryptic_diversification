################################################################################
#### PHYLOGENETIC INFERENCE ####
################################################################################

scripts_dir=/home/nibtve93/scripts/phylogeneticInference

prefix=MicrocebusPhylogenomics
set_id=fullSet
vcf_dir=$PWORK/$prefix/gatk
phyl_dir=$PWORK/$prefix/phylogeneticInference

mkdir -p $phyl_dir/logFiles

#################################################################
#### 1 MAXIMUM LIKELIHOOD INFERENCE WITH ASCERTAINMENT BIAS CORRECTION ####
#################################################################
mkdir -p $phyl_dir/ml

## Convert VCF file to PHYLIP format and calculate basic alignment statistics for different maximum missing data thresholds (ranging from 0.05% to 0.95%)
format=phylip # Output format of alignment (phylip or nexus)
for i in 0.05 0.25 0.5 0.75 0.95
do
format=phylip # Output format of alignment (phylip or nexus)
sbatch --account=nib00015 --output=$phyl_dir/logFiles/vcf_convert.ml.$format.$set_id.maxmiss$i.oe $scripts_dir/vcf_convert.sh $scripts_dir $vcf_dir/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$i.vcf $phyl_dir/ml $format
done

## Run phylogenetic inference with ascertainment bias correction in IQ-TREE
njobs=12 # Since there was walltime limit of 48 h on the server, we used DMTCP (https://dmtcp.sourceforge.io/) for checkpointing; $njobs jobs were submitted, each of which would only run a maximum of 48 h 
ufboot=1000

for i in 0.05 0.25 0.5 0.75 0.95
do
out_dir=$phyl_dir/ml/maxmiss$i
mkdir -p $out_dir

for j in $(seq 1 $njobs)
do
echo -e "#### Submitting job $i of missing data threshold $i"
# If first script:
if [ $j == 1 ]
then
# Declare directory to save checkpoint files (will be created in submitted script)
check_out=$out_dir/checkpoint_$j
# Submit job and save submission ID
jid=$(sbatch --account=nib00015 --output=$out_dir/iqtree.$set_id.maxmiss$i.job$j.oe $scripts_dir/iqtree_checkpoint.sh $phyl_dir/ml/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$i.noinv.phy $out_dir/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$i.noinv $ufboot $check_out)
declare runid_$j=${jid##* }

# If not first script:
else
# Assign input directory (which is output directory of previous iteration)
check_in=$check_out
# Declare directory to save checkpoint files
check_out=$out_dir/checkpoint_$j
# Get submission ID of previous iteration
varname=runid_$(( $j - 1 ))
# Submit next job and save submission ID
jid=$(sbatch --account=nib00015 --output=$out_dir/iqtree.$set_id.maxmiss$i.job$j.oe --dependency=afterany:${!varname} $scripts_dir/iqtree_checkpoint_cont.sh $check_in $check_out)
#see above
declare runid_$j=${jid##* }
fi
done

done

#################################################################
#### 2 QUARTET-BASED INFERENCE FOR INDIVIDUAL AND POPULATION ASSIGNMENT ####
#################################################################
mkdir -p $phyl_dir/quartet

for i in 0.05 0.25 0.5 0.75 0.95
do
## Thin VCF file by 10,000 bp intervals to ensure independence of SNPs
vcftools --vcf $vcf_dir/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$i.vcf --thin 10000 --recode --recode-INFO-all --stdout > $vcf_dir/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$i.thin10k.vcf

## Convert VCF file to NEXUS format and calculate basic alignment statistics
format=nexus # Output format of alignment (phylip or nexus)
sbatch --account=nib00015 --output=$phyl_dir/logFiles/vcf_convert.quartet.$format.$set_id.maxmiss$i.oe $scripts_dir/vcf_convert.sh $scripts_dir $vcf_dir/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$i.thin10k.vcf $phyl_dir/quartet $format

## Create taxon partitions block files
# For population assignment, the file was created manually and saved as $phyl_dir/quartet/$set_id.maxmiss$i.paup.population.nex
# For individual assignment:
echo "BEGIN SETS;" > $phyl_dir/quartet/$set_id.maxmiss$i.taxPartitions.individual.nex
echo -e "\t TAXPARTITION SPECIES =" >> echo "BEGIN SETS;" > $phyl_dir/quartet/$set_id.maxmiss$i.taxPartitions.individual.nex
for j in $(awk '{ print $1 }' $phyl_dir/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$i.noinv.phy | tail -n+2)
do
	echo -e "\t\t${j}:${j}," >> echo "BEGIN SETS;" > $phyl_dir/quartet/$set_id.maxmiss$i.taxPartitions.individual.nex
done
sed -i '$ s/,$//' echo "BEGIN SETS;" > $phyl_dir/quartet/$set_id.maxmiss$i.taxPartitions.individual.nex # Remove last comma
echo -e "\t\t;" >> echo "BEGIN SETS;" > $phyl_dir/quartet/$set_id.maxmiss$i.taxPartitions.individual.nex
echo "END;" >> echo "BEGIN SETS;" > $phyl_dir/quartet/$set_id.maxmiss$i.taxPartitions.individual.nex
echo "" >> echo "BEGIN SETS;" > $phyl_dir/quartet/$set_id.maxmiss$i.taxPartitions.individual.nex

## Create PAUP block files
nt=80
nq=20000000
seed=$RANDOM
# For population assignment:
echo "BEGIN PAUP;" > $phyl_dir/quartet/$set_id.maxmiss$i.paup.population.nex
echo -e "\toutgroup SPECIES.Cheirogaleusmaj_106245 SPECIES.Cheirogaleusmed_106354 SPECIES.Cheirogaleuscro_106189;" >> $phyl_dir/quartet/$set_id.maxmiss$i.paup.population.nex
echo -e "\tset root=outgroup outroot=monophyl;" >> $phyl_dir/quartet/$set_id.maxmiss$i.paup.population.nex
echo -e "\tsvdq nthreads=$nt nquartets=${nq} evalQuartets=random taxpartition=SPECIES bootstrap=standard seed=$seed;" >> $phyl_dir/quartet/$set_id.maxmiss$i.paup.population.nex
echo -e "\tsavetrees format=Newick file=$phyl_dir/quartet/final.maxmiss$i.thinned.speciesAssignment.svdq.tre savebootp=nodelabels;" >> $phyl_dir/quartet/$set_id.maxmiss$i.paup.population.nex
echo -e "\tquit;" >> $phyl_dir/quartet/$set_id.maxmiss$i.paup.population.nex
echo "END;" >> $phyl_dir/quartet/$set_id.maxmiss$i.paup.population.nex
# For individual assignment
echo "BEGIN PAUP;" > $phyl_dir/quartet/$set_id.maxmiss$i.paup.individual.nex
echo -e "\toutgroup SPECIES.Cheirogaleusmaj_106245 SPECIES.Cheirogaleusmed_106354 SPECIES.Cheirogaleuscro_106189;" >> $phyl_dir/quartet/$set_id.maxmiss$i.paup.individual.nex
echo -e "\tset root=outgroup outroot=monophyl;" >> $phyl_dir/quartet/$set_id.maxmiss$i.paup.individual.nex
echo -e "\tsvdq nthreads=$nt nquartets=${nq} evalQuartets=random taxpartition=SPECIES bootstrap=standard seed=$seed;" >> $phyl_dir/quartet/$set_id.maxmiss$i.paup.individual.nex
echo -e "\tsavetrees format=Newick file=$phyl_dir/quartet/final.maxmiss$i.thinned.svdq.tre savebootp=nodelabels;" >> $phyl_dir/quartet/$set_id.maxmiss$i.paup.individual.nex
echo -e "\tquit;" >> $phyl_dir/quartet/$set_id.maxmiss$i.paup.individual.nex
echo "END;" >> $phyl_dir/quartet/$set_id.maxmiss$i.paup.individual.nex

## Concatenate files and submit SVDquartets job
for j in population individual
do
	cat $phyl_dir/quartet/allScaffolds.annot.SNP.minInd.DP.mac.GATKfilt-hard.maxmiss$i.thin10k.noinv.nex $phyl_dir/quartet/$set_id.maxmiss$i.taxPartitions.$j.nex $phyl_dir/quartet/$set_id.maxmiss$i.paup.$j.nex > $phyl_dir/quartet/$set_id.maxmiss$i.paup.$j.concat.nex
	sbatch --output=$phyl_dir/logFiles/svdq.$set_id.maxmiss$i.oe $scripts_dir/svdq.sh $phyl_dir/quartet/$set_id.maxmiss$i.paup.$j.concat.nex $phyl_dir/quartet/$set_id.maxmiss$i.paup.$j.concat.nex.log
done

done
