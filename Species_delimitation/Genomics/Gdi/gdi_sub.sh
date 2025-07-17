################################################################################
#### SPECIES DELIMITATION THROUGH GENEALOGICAL DIVERGENCE (GDI) ####
################################################################################

scripts_dir=/home/nibtve93/scripts/gdi

prefix=MicrocebusPhylogenomics

div_dir=$PWORK/$prefix/divergenceTimeEstimation
gdi_dir=$PWORK/$prefix/speciesDelimitation/genomics/gdi
lin_map=$div_dir/lineager_map.txt # file that maps node ids to lineages (can be inferred from bpp --summary output)

#################################################################
#### INFERENCE OF GENEALOGICAL DIVERGENCE INDEX ####
#################################################################
conda activate R

burnin=100000 # Burn-in
last_sample=900000 # Last MCMC sample to be considered
mutrate=1.236 # Gamma distribution of mutation rate will have mean $mutrate * 10e-8
mutrate_var=0.107 # Gamma distribution of mutation rate will have variance $mutrate_var * 10e-8
gentime=3.5 # Lognormal distribution of generation time will have mean ln($gentime)
gentime_sd=1.16 # Lognormal distriubtion of generation time will have standard deviation ln($gentime_sd)

## Convert MCMCM files
replicates=4

for j in $(seq 1 $replicates)
do
echo "... Processing run $j"
mkdir -p $div_dir/run_$j
while read -r id lineage; do
  sed -e "s/theta_${id}/theta_${lineage}/g" -e "s/tau_${id}/tau_${lineage}/g" $div_dir/run_$j/mcmc.run$j.txt > $gdi_dir/run_$j/mcmc.run$j.txt
done < "$lin_map"
done

## Estimate genealogical divergence index (gdi)
for j in $(seq 1 $replicates)
do
echo "... Processing run $j"
Rscript $scripts_dir/process_logs.R $gdi_dir/run_$j/mcmc.run$j.txt $gdi_dir/run_$j/gdi.run$j.txt $burnin $last_sample
done

## Average gdi estimates
mkdir -p $gdi_dir/run_average
chains=()
for j in $(seq 1 $replicates)
do
  chains+=("$gdi_dir/run_$j/gdi.run$j.txt")
done
Rscript $scripts_dir/average_chains.R "${chains[@]}" $gdi_dir/run_average/gdi.average.txt
