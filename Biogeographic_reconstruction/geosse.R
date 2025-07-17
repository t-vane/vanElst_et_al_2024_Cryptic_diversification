#!/usr/bin/env Rscript

### Script adapted and modified from Everson et al. (2016), Syst. Biol.

# Libraries
library(diversitree)
library(phytools)
library(geiger)

# Load ultrametric tree
tree <- read.tree("tree.nwk")
drop_tips <- c("Mmittermeieri", "MlehilahytsaraS", "Msp1", "Mboraha", "Mmarohita", "Mmanitatra", "Mganzhorni", "MmurinusC", "Mbongolavensis")
tree <- drop.tip(tree, tip = drop_tips)

# Load geo data to tip states on phylogeny#
tree$tip.state<-c(
  "Mdanfossi"=1,
  "Mravelobensis"=1,
  "Mgriseorufus"=1,
  "MmurinusN"=1,
  "Mjonahi"=2,
  "Mmacarthurii"=2,
  "Mgerpi"=2,
  "Mjollyae"=2,
  "Msimmonsi"=2,
  "Mtavaratra"=1,
  "Marnholdi"=2,
  "Msambiranensis"=2,
  "Mmargotmarshae"=2,
  "Mmamiratra"=2,
  "Mtanosi"=2,
  "MlehilahytsaraN"=2,
  "Mmyoxinus"=1,
  "Mberthae"=1,
  "Mrufus"=2
  )
statecols=c("violet","#FFD27F","#7FB17E")

#######################################
####### MAXIMUM LIKELIHOOD ############
#######################################

# Get starting parameters for GeoSSE analysis
p<-starting.point.geosse(tree)

## Model construction
# Full model
full<-make.geosse(tree, tree$tip.state, strict=FALSE)
# Model without speciation in state 0
nosAB<-constrain(full, sAB ~ 0)
# Model with no extinction
noEx<-constrain(full, xA ~ 0, xB ~ 0)
# Model with no dispersal
noD<-constrain(full, dA ~ 0, dB ~ 0)
# Model without regional dependence of speciation
sameS<-constrain(full, sA ~ sB)
# Model with no regional difference in dispersal
sameD<-constrain(full, dA ~ dB)
# Model with no regional difference in extinction
sameEx<-constrain(full, xA ~ xB)
# Model with no extinction and no regional difference in speciation
noExsameS<-constrain(full, sA ~ sB, xA ~ 0, xB ~ 0)
# Model with no extinction and same dispersal
noExsameD<-constrain(full, xA ~ 0, xB ~ 0, dA ~ dB)
# Model with no extinction, same dispersal, same speciation
noExsameSsameD<-constrain(full, xA ~ 0, xB ~ 0, sA ~ sB, dA ~ dB)

## ML parameter estimation for each model
# Full model
print("... Estimating full likelihood model")
MLfull<-find.mle(full, p)
save.image(file = "Mlfull.RData")
p<-coef(MLfull)

# Define constrained models
models <- c("nosAB", "noEx", "noD", "sameS", "sameD", "sameEx", "noExsameS", "noExsameD", "noExsameSsameD")
# Loop through constrained models
for (model_name in models) {
  print(paste("... Estimating", model_name, "model"))
  
  model_obj <- get(model_name)
  mle_result <- find.mle(model_obj, p[argnames(model_obj)])
  
  assign(paste0("ML", model_name), mle_result)
  save.image(file = paste0(model_name, ".RData"))
}

## Model comparison
print("... Comparing models")
round(rbind(full=coef(MLfull), 
            no.sAB=coef(MLnosAB, TRUE), 
            no.ex=coef(MLnoEx, TRUE), 
            same.S=coef(MLsameS, TRUE), 
            same.D=coef(MLsameD, TRUE), 
            same.Ex=coef(MLsameEx, TRUE),
            noEx.sameS=coef(MLnoExsameS, TRUE), 
            noEx.sameD=coef(MLnoExsameD, TRUE), 
            noEx.sameS.sameD=coef(MLnoExsameSsameD, TRUE)), 
            9)
anova(MLfull, 
      no.sAB=MLnosAB, 
      no.ex=MLnoEx, 
      same.S=MLsameS, 
      same.D=MLsameD, 
      same.Ex=MLsameEx,
      noEx.sameS=MLnoExsameS, 
      noEx.sameD=MLnoExsameD, 
      noEx.sameS.sameD=MLnoExsameSsameD)

AICvals<-c(267.95,
           269.76,
           269.89,
           271.61,
           271.89,
           271.95,
           273.89,
           275.38,
           281.06)
akaike.weights(AICvals)
save.image(file = "final.RData")
