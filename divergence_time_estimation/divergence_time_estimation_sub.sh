################################################################################
#### DIVERGENCE TIME ESTIMATION ####
################################################################################

scripts_dir=/home/nibtve93/scripts/divergenceTimeEstimation

prefix=MicrocebusPhylogenomics
set_id=reducedSet

fasta_dir=$PWORK/$prefix/gatk/locusExtraction/$set_id/fasta/"$set_id"_bylocus_final
div_dir=$PWORK/$prefix/divergenceTimeEstimation
subset_dir=$div_dir/subset
seq_file_all=$div_dir/seq_file_all.txt
seq_file_6k=$div_dir/seq_file_6k.txt
locus_dir=$div_dir/loci

mkdir -p $div_dir/logFiles

#################################################################
#### 1 PREPARE FULL BPP SEQUENCE FILE FOR BPP ####
#################################################################
mkdir -p $locus_dir

## Create per-locus files in correct format
sbatch --account=nib00015 --output=$div_dir/logFiles/prepare_loci.oe $scripts_dir/prepare_loci.sh $fasta_dir $locus_dir $prefix

## Concatenate to monolithic file
find $locus_dir -maxdepth 1 -name "*locus" -type f -exec cat {} + >> $seq_file_all

#################################################################
#### 2 PREPARE SUBSAMPLED BPP SEQUENCE FILE (6000 LOCI WITH LEAST AMOUNT OF MISSING DATA) ####
#################################################################
mkdir -p $subset_dir/seqFiles/

## Create per-locus seq files
n_loci=$(( $(cat $seq_file_all | wc -l) / 58 ))
for i in $(seq 1 $n_loci)
do
j=$(( $i - 1 ))
echo $i
sed -n "$(( 1 + $j * 58 )),$(( 58 + $j * 58))p" $seq_file_all > $subset_dir/seqFiles/seq_file_locus_$i.txt
done 

## Calculate missingness for each locus
rm $subset_dir/loci.missing.stats
for i in $(seq 1 $n_loci)
do
echo $i
n_all=$(awk '{print $2}' $subset_dir/seqFiles/seq_file_locus_$i.txt | tail -n +2 | wc -c)
n_miss=$(awk '{print $2}' $subset_dir/seqFiles/seq_file_locus_$i.txt | tail -n +2 | awk -v RS='' -v FPAT='N' '{print NF}')
p_miss=$(lua -e "print($n_miss/$n_all)")
echo -e "locus_$i\t$p_miss" >> $subset_dir/loci.missing.stats
done

## Extract the 6000 locus ids and create seq file
sort -k 2 $subset_dir/loci.missing.stats > $subset_dir/loci.missing.sorted.stats | head -n 6000 | awk '{print $1}' > $subset_dir/loci_6k.txt

rm $div_dir/seq_file.6k.txt
for i in $(cat $subset_dir/loci_6k.txt)
do
echo $i
cat $subset_dir/seqFiles/seq_file_locus_$i.txt >> $seq_file_6k
done

#################################################################
#### 3 CONFIRM TOPOLOGY OF SELECTED LOCI SUPPORTS SPECIES TREE ####
#################################################################
mkdir -p $subset_dir/phylogeneticInference/fasta

## Convert to FASTA
for i in $(cat $subset_dir/loci_6k.txt)
do
echo $i
tail -n +2 $subset_dir/seqFiles/seq_file_locus_$i.txt | sed -E 's/\^/>/g; s/\s+/\n/g' > $subset_dir/phylogeneticInference/fasta/locus_$i.txt
done

## Concatenate loci
AMAS.py concat -i $subset_dir/phylogeneticInference/fasta/locus_*.txt -f fasta -d dna -u phylip -t $subset_dir/phylogeneticInference/alignment.concat.phy

## Run IQTREE
sbatch --output=$div_dir/logFiles/iqtree_bpp.oe $scripts_dir/iqtree_bpp.sh $subset_dir/phylogeneticInference/alignment.concat.phy $subset_dir/phylogeneticInference/alignment.concat 1000

#################################################################
#### 4 RUN BPP ####
#################################################################
NJOBS=15 # number of checkpointing runs per independent bpp run
replicates=4

for j in $(seq 1 $replicates)
do
out_dir=$div_dir/run_$j
ctrl_file=$out_dir/A00.bpp.6k.40threads.run$j.ctl

mkdir -p $out_dir

for i in $(seq 1 $njobs)
do
echo -e "... Submitting job $i of run $j"
# If first script:
if [ $i == 1 ]
then
# Declare directory to save checkpoint files
check_out=$out_dir/checkpoint_$i
mkdir -p $check_out
# Submit job and save submission ID
jid=$(sbatch --account=nib00015 --output=$out_dir/$i.bpp.oe $scripts_dir/bpp_checkpoint.sh $ctrl_file $check_out)
declare runid_$i=${jid##* }

# If not first script:
else
# Assign input directory (which is output directory of previous iteration)
check_in=$check_out
# Declare directory to save checkpoint files
check_out=$out_dir/checkpoint_$i
# Get submission ID of previous iteration
varname=runid_$(( $i - 1 ))
# Submit next job and save submission ID
jid=$(sbatch --account=nib00015 --output=$out_dir/$i.bpp.oe --dependency=afterany:${!varname} $scripts_dir/bpp_checkpoint_rerun.sh $check_in $check_out)
declare runid_$i=${jid##* }
fi
done

done

#################################################################
#### 5 POST-PROCESSING ####
#################################################################
conda activate R

burnin=100000 # Burn-in
last_sample=900000 # Last MCMC sample to be considered
mutrate=1.236 # Gamma distribution of mutation rate will have mean $mutrate * 10e-8
mutrate_var=0.107 # Gamma distribution of mutation rate will have variance $mutrate_var * 10e-8
gentime=3.5 # Lognormal distribution of generation time will have mean ln($gentime)
gentime_sd=1.16 # Lognormal distriubtion of generation time will have standard deviation ln($gentime_sd)

## Convert MCMC estimates of single chains to divergence time in years and effective population size and create summary trees
for j in $(seq 1 $replicates)
do
echo "... Processing run $i"
mcmc_file=$div_dir/run_$j/mcmc.run$i.txt

for i in dist nodist
do
[ "$i" = "dist" ] && dist=TRUE || dist=FALSE
Rscript $scripts_dir/convert_mcmc.R $div_dir/run_$j $mcmc_file $div_dir/run_$j/mcmc.run$j.conv.$i.txt $mutrate_gen $mutrate_var $gentime $gentime_sd $burnin $last_sample $dist
done

# Create tree
bpp --summary $out_dir/A00.bpp.6k.40threads.run$j.ctl
done


## Average chains and create tree
for i in dist nodist
do
# Average chains
mkdir -p $div_dir/run_average
chains=()
for j in $(seq 1 $replicates)
do
  chains+=("$div_dir/run_$j/mcmc.run$j.conv.$i.txt")
done
Rscript $scripts_dir/average_chains.R "${chains[@]}" $div_dir/run_average/mcmc.mean.convert.$i.txt

# Create tree
bpp --summary $out_dir/A00.bpp.6k.40threads.mean.$i.ctl
done
