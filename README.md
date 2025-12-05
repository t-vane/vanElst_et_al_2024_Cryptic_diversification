# Integrative taxonomy clarifies the evolution of a cryptic primate clade

This repository holds scripts for the following analyses conducted as part of the publication [van Elst et al. (2025), *Nat. Ecol. Evol.*](https://doi.org/10.1038/s41559-024-02547-w):
- Genotyping
- Phylogenetic inference
- Species delimitation
- Divergence time estimation
- Biogeographic reconstruction
- Modelling morphological and climatic niche evolution
- Plots

Input and output files can be found in the [Dryad digital repository](https://doi.org/10.5061/dryad.b2rbnzsp3). 

### Genotyping
`./genotyping` contains the following subdirectories:

`./reference_mapping` contains scripts to align cleaned reads against the *Microcebus murinus* reference genome (Mmur 3.0; [Larsen et al. (2017), *BMC Biol.*](https://doi.org/10.1186/s12915-017-0439-6) with [BWA v0.7.17-r1198-dirty](https://github.com/lh3/bwa) and to filter BAM files with [SAMtools v1.11](http://www.htslib.org/).

`./genotype_calling` contains scripts to call genotypes with [GATK v4.1.9.0](https://gatk.broadinstitute.org/hc/en-us) and apply various filters to VCF files with [VCFtools v0.1.17](https://vcftools.github.io/index.html), [Bcftools v1.11](https://samtools.github.io/bcftools/) and [GATK v3.8.1/v4.1.9.0](https://gatk.broadinstitute.org/hc/en-us).

`./locus_extraction` contains scripts to convert called genotypes to phased RAD loci following [Poelstra et al. (2021), *Syst. Biol.*](https://doi.org/10.1093/sysbio/syaa053)

`./Genotype_likelihoods` contains scripts to infer genotype likelihoods for candidate species subsets with [ANGSD v0.92](http://www.popgen.dk/angsd/index.php/ANGSD).

### Phylogenetic inference
`./phylogenetic_inference` contains scripts to infer maximum likelihood and quartet-based phylogenies with [IQ-TREE v2.2.0](http://www.iqtree.org/#download) and [SVDquartets of PAUP* v4.0a](https://paup.phylosolutions.com/), respectively.

### Species delimitation
`./species_delimitation` contains the following subdirectories:

`./genomics` contains scripts to infer individual ancestries with [NGSadmix v32](https://www.popgen.dk/software/index.php/NgsAdmix) and to estimate of genealogical divergence indices (*gdi*) between candidate species. Scripts to test for genomic patterns of isolation-by-distance can be found [here](https://github.com/gsgarlata/Isolation-by-distance-within-vs-between-species/tree/main).

`./morphometry` contains scripts to estimate overlap in morphometry in the R package ['dynRB' v0.18](https://cran.r-project.org/web/packages/dynRB/index.html) and test for morphometric patterns of isolation-by-distance.

`./climatic_niche` contains scripts to infer climatic niche models and estimate differentiation among these with the R package ['ENMtools' v1.0.7](https://github.com/danlwarren/ENMTools).

`./acoustic_communication` contains scripts to estimate differentiation in acoustic communication.

`./reproductive_activity` contains scripts to estimate differentiation in reproductive activity.

### Divergence time estimation
`./divergence_time_estimation` contains scripts to infer divergence times from phased RAD loci under a coalescent model in [BPP v4.4.1](https://github.com/bpp/bpp).

### Biogeographic reconstruction
`./biogeographic_reconstruction` contains scripts to infer ancestral biogeogeography in the R package ['BioGeoBEARS' v1.1.2](https://github.com/nmatzke/BioGeoBEARS) and habitat-associated speciation rates in the R package ['diversitree' v0.9-16](https://www.zoology.ubc.ca/prog/diversitree/).

### Modelling morphological and climatic niche evolution
`./modelling_morphological_and_climatic_niche_evolution` contains the following subdirectories:

`./climatic_niche_evolution` contains scripts to model the evolution of climatic niches along the inferred *Microcebus* phylogeny.

`./morphological_evolution` contains scripts to model the evolution of morphology along the inferred *Microcebus* phylogeny, including cross-validation.

Modelling was conducted with the R packages ['phytools' v2.3-0](https://cran.r-project.org/web/packages/phytools/index.html), ['dynRB' v0.18](https://cran.r-project.org/web/packages/dynRB/index.html), ['phyloclim' v0.9.5](https://cran.r-project.org/web/packages/phyloclim/index.html), ['mvMORPH' v1.1.9](https://github.com/JClavel/mvMORPH), ['tmvtnorm' v1.6](https://github.com/cran/tmvtnorm).

### Plots
`./plots` contains scripts for various plots that were generated in R for the publication.