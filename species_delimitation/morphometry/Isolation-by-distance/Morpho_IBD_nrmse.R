### Script by Jordi Salmona

##### load libraries ####
library(sp)
library(usedist)
library(ggplot2)
library(ggrepel)
library(vegan)
library(plyr)
library(dplyr)
library(scales)
library(ggplot2)
library(reshape2)
library(readxl)         # data import
library(xlsx)
library(MASS)           # LDA
library(ade4)           # Mantel.randtest
library(viridis)
library(ggridges)
library(ggpubr)

#### set working directory #####
setwd("D:/jordi/mada/consortium_rad_microcebus/js_phylo/IBD/IBD_morpho")
data <- read_excel("D:/jordi/mada/consortium_rad_microcebus/js_phylo/IBD/IBD_morpho/data_step5.xlsx", sheet = "data", col_names = TRUE)

#### set color codes #### 
# depends on whether you compare 3 taxa or 2 taxa
colorsFor3 <- c("#0072B2", "#D55E00", "#CC79A7", "#E69F00","#56B4E9","#009E73")
shapesFor3 <- c(3,4,8,15,16,17)
colorsFor2 <- c("#CC79A7", "#E69F00", "#56B4E9")
shapesFor2 <- c(8,15,16)

#### data preparation ####
str(data)
head(data)
data$population <- as.factor(data$population)

# rename to fit 4 letter names for later use
levels(data$population)
data$population <- recode(data$population, gris_Beza = 'gris', 
                          lehi_S = 'lehi',
                          lehi_N = 'lehi',
                          lehi_mitt = 'mitt',
                          sp1 = "Msp1",
                          muri_C = "murC",
                          muri_N = "murN",
                          maro = "joll")

#### define species list ####
species_list = list(c('mitt','lehi'),c('tavN','tavC'),c('tavC','tavS'),c('tavN','tavS'),
                    c('rave','danf'),c('rave','bong'),c('bong','danf'),
                    c('simm','bor'),
                    c('mami','marg'), c('mami','samb'),c('marg','samb'),
                    c('myox','bert'), c('bert','rufu'),c('myox','rufu'), #rufus has only 5 variables I need to twist the code
                    c('ganz','mani'),c('ganz','murC'),c('ganz','murN'), c('murC','murN'),
                    c('maca','jona'),
                    c('gerp','joll'),#c('maro','gerp'),c('joll','maro'),
                    c('arno','Msp1'))
species_list

#### old stuff ####
# extract all taxa to single dataframes, to look if downsampling is needed; combine after; e.g., bora, maro not possible
bong <- dplyr::filter(data, population == "bong")
danf <- dplyr::filter(data, population == "danf")
ganz <- dplyr::filter(data, population == "ganz")
gerp <- dplyr::filter(data, population == "gerp")
gris <- dplyr::filter(data, population == "gris")
joll <- dplyr::filter(data, population == "joll")
jona <- dplyr::filter(data, population == "jona")
mitt <- dplyr::filter(data, population == "mitt")
lehi <- dplyr::filter(data, population == "lehi")
maca <- dplyr::filter(data, population == "maca")
mami <- dplyr::filter(data, population == "mami")
marg <- dplyr::filter(data, population == "marg")
murC <- dplyr::filter(data, population == "murC")
murN <- dplyr::filter(data, population == "murN")
myox <- dplyr::filter(data, population == "myox")
rave <- dplyr::filter(data, population == "rave")
rufu <- dplyr::filter(data, population == "rufu")
samb <- dplyr::filter(data, population == "samb")
simm <- dplyr::filter(data, population == "simm")
Msp1 <- dplyr::filter(data, population == "Msp1")
tano <- dplyr::filter(data, population == "tano")
tava <- dplyr::filter(data, population == "tava")

#murN <- sample_n(murN, 150) # re-sampling to reduce sample size
#rave <- sample_n(rave, 100) 
#tava <- sample_n(tava, 100) 


#### start of the species pairs' loop ####
#spnb<-1
species_list<-list() ; for(spnb in 1:length(ibd_all_dat[,1])){ species_list[[spnb]]<-c(ibd_all_dat[spnb,1],ibd_all_dat[spnb,2])}
species_list

for(spnb in 17:length(species_list)){
#for(spnb in c(20)){ #:14,16
print(species_list[[spnb]])
sp1<-species_list[[spnb]][[1]]
sp2<-species_list[[spnb]][[2]]

  #for(spnb in 3:4){
  print(paste0("processing pair ",spnb))
  
  
#### define duo/triplet for analysis ####
dat_sp1 <- dplyr::filter(data, population == sp1)
dat_sp2 <- dplyr::filter(data, population == sp2)

if(dim(dat_sp1)[1]<5 | dim(dat_sp2)[1]<5){print(paste("not enough data for pair",species_list[[spnb]][[1]],"-",
                                                      species_list[[spnb]][[2]],"moving to next pair"))}
if(dim(dat_sp1)[1]>=5 & dim(dat_sp2)[1]>=5){print(paste("enough data for pair",species_list[[spnb]][[1]],"-",
                                                      species_list[[spnb]][[2]],"processing"))
data1 <- rbind(dat_sp1,dat_sp2) 

#### select morphological variables ####
# select variables 7 var with significant differences between taxa in ANOVA
#data2 <- dplyr::select(data1, population, ear.length, head.length, head.width, interorbital.dist,snout.length, lower.leg.length, hind.foot.length,latitude, longitude)
data2 <- data1[,c("population", "ear.length", "head.length", "head.width", 
                  "interorbital.dist","snout.length", "lower.leg.length", 
                  "hind.foot.length","latitude", "longitude")]

# select variables 12 (all) variables
variables<-c("population", "ear.length", "ear.width",	"head.length",	"head.width",
             "interorbital.dist",	"intraorbital.dist",	"snout.length",	"lower.leg.length",	"hind.foot.length",
             "third.toe.length",	"body.length",	"tail.length","body.mass","latitude", "longitude")
data2_12var <- data1[,variables]
data2_12var_sp1<- dat_sp1[,variables]
data2_12var_sp2<- dat_sp2[,variables]

dim(data2_12var)
print("which variables are fully missing in each of species candidates ?")
mdat<-c(which((colSums(is.na(data2_12var_sp1))==dim(data2_12var_sp1)[1])==T),
            which((colSums(is.na(data2_12var_sp2))==dim(data2_12var_sp2)[1])==T))
print(mdat)
if(length(mdat)>=1){print(paste(length(mdat),"variables is/are missing in at least 1 sp and removed for both"))
  data2_12var <- data2_12var[,which(colnames(data2_12var)!=unique(names(mdat)))]
  print(dim(data2_12var))
  data2 <- data2[,which(colnames(data2)!=unique(names(mdat)))]
}

print(colSums(is.na(data2)))
print(plyr::count(colSums(is.na(data2))))
print(plyr::count(rowSums(is.na(data2))))
print(colSums(is.na(data2_12var)))
print(plyr::count(colSums(is.na(data2_12var))))
print(plyr::count(rowSums(is.na(data2_12var))))

#### select missingness trhesholds ####
# discard all missing data
#data4 <- na.omit(data2)
# keep all missing data
data4_na <- droplevels(data2)
data4_12var <- droplevels(data2_12var)
# keep only samples with max 30% of missing data 
nvar<-(which(colnames(data2)=="latitude")-2)
data4_na30<-data2[rowSums(is.na(data2))<(0.3*nvar),]
nvar12<-(which(colnames(data2_12var)=="latitude")-2)
data4_12var_na30 <- data2_12var[rowSums(is.na(data2_12var))<(0.3*nvar12),]


#### verify data length before lauching loop ####
if(dim(data4_na30[data4_na30[,1]==sp1,])[1]<5 | dim(data4_na30[data4_na30[,1]==sp2,])[1]<5){
  print(paste("not enough data for pair",species_list[[spnb]][[1]],"-",species_list[[spnb]][[2]],"moving to next pair"))
  } 
if(dim(data4_12var_na30[data4_12var_na30[,1]==sp1,])[1]<5 | dim(data4_12var_na30[data4_12var_na30[,1]==sp2,])[1]<5){
  print(paste("not enough data for pair",species_list[[spnb]][[1]],"-",species_list[[spnb]][[2]],"moving to next pair"))}


data_list<-list(data4_na30,data4_12var_na30)
data_types<-c("7var_30NA","12var_30NA")
#head(data_list[[1]])
last_var<-c((which(colnames(data_list[[1]])=="latitude")-1),(which(colnames(data_list[[2]])=="latitude")-1))
if(dim(data4_na30[data4_na30[,1]==sp1,])[1]<5 | dim(data4_na30[data4_na30[,1]==sp2,])[1]<5){
  print("not enough data with 7 variables keeping only the larger dataset")
  data_list<-list(data4_12var_na30)
  data_types<-c("12var_30NA")
  #head(data_list[[1]])
  last_var<-c((which(colnames(data_list[[1]])=="latitude")-1))
} 
if(dim(data4_12var_na30[data4_12var_na30[,1]==sp1,])[1]<=5 | dim(data4_12var_na30[data4_12var_na30[,1]==sp2,])[1]<=5){
  print("not enough data with the larger dataset building a default dataset with available data")
  print("which variables are available in sp1 and then in sp2")
  print(colSums(is.na(data2_12var[data2_12var[,1]==sp1,])))
  print(colSums(is.na(data2_12var[data2_12var[,1]==sp2,])))
  # keep only variables present in at least 30% of the samples of each species 
  nsamp1<-dim(data2_12var[data2_12var[,1]==sp1,])[1] ; print(nsamp1)
  nsamp2<-dim(data2_12var[data2_12var[,1]==sp2,])[1] ; print(nsamp2)
  rm(data4_12var_na30)
  data4_12var_na30 <- data2_12var[,colSums(is.na(data2_12var[data2_12var[,1]==sp1,]))<(0.3*nsamp1) & 
                colSums(is.na(data2_12var[data2_12var[,1]==sp2,]))<(0.3*nsamp2)]
  nvar12<-(which(colnames(data4_12var_na30)=="latitude")-2)
  data4_12var_na30 <- data4_12var_na30[rowSums(is.na(data4_12var_na30))<(0.3*nvar12),]
  data4_12var_na30 <- droplevels(data4_12var_na30)
  data_list<-list(data4_12var_na30)
  data_types<-c("12var_30NA")
  #head(data_list[[1]])
  last_var<-c((which(colnames(data_list[[1]])=="latitude")-1))
  
}
if(dim(data4_12var_na30[data4_12var_na30[,1]==sp1,])[1]>=5 & dim(data4_12var_na30[data4_12var_na30[,1]==sp2,])[1]>=5){

#### add populations and ids ####
for(i in 1:length(data_list)){ 
  data_list[[i]]<-droplevels(data_list[[i]])
  data_list[[i]][,2:last_var[[i]]]<-scale(data_list[[i]][,2:last_var[[i]]])
  data_list[[i]]$consecutive_numbers <- 100 + 1:nrow(data_list[[i]]) 
  data_list[[i]]$spp_number <- paste(data_list[[i]]$population, data_list[[i]]$consecutive_numbers, sep="_") 
  data_list[[i]]$latitude <- as.numeric(data_list[[i]]$latitude)
  data_list[[i]]$longitude <- as.numeric(data_list[[i]]$longitude)
  data_list[[i]]$spp_number <- as.factor(data_list[[i]]$spp_number)
  print(table(data_list[[i]]$population))}  
  
#### generate morphometric distance matrices ####
mat_list<-list()
distm=c("euclidean","maximum","manhattan","canberra")
for(i in 1:length(data_list)){
  mat_list[[i]]<-list()
  for(j in 1:length(distm)){
    mat_list[[i]][[j]]<-dist(data_list[[i]][,2:last_var[[i]]],method=distm[j])
    mat_list[[i]][[j]] <- dist_setNames(mat_list[[i]][[j]], data_list[[i]]$spp_number)}}
save(mat_list,file=paste0(sp1,"_",sp2,"_matlist.rda"))    

#### generate geographic distance matrices ####
geo_list<-list()
for(i in 1:length(data_list)){ 
  geo_list[[i]]<-as.dist(spDists(as.matrix(data_list[[i]][,c("longitude","latitude")]),longlat = T))
  geo_list[[i]][which(geo_list[[i]]==0)]<-rnorm(length(which(geo_list[[i]]==0)),1,0.3)}
save(geo_list,file=paste0(sp1,"_",sp2,"_geolist.rda"))

#### turn dist matrices to tables with ind-pairs, species and categories ####
print("turn dist matrices to tables with ind-pairs, species and categories")
df_list<-list()
for(i in 1:length(data_list)){
  df_list[[i]]<-list()
  for(j in 1:length(distm)){
    indPairs<-c()
    morphoD <- c()
    geoD <- c()
    morphomat<-as.matrix(mat_list[[i]][[j]])
    geomat<-as.matrix(geo_list[[i]])
    for(row in 1:nrow(morphomat)){
      for(col in 1:ncol(morphomat)){
        if (row < col) {
          indPairs <- c(indPairs,paste0(rownames(morphomat)[row],"_",
                                        colnames(morphomat)[col]))
          morphoD <- c(morphoD, morphomat[row, col])
          geoD <- c(geoD, geomat[row, col])
    }}}
    combined <- data.frame(indPairs, morphoD, geoD)
    combined$species1 <-substr(combined[,1], 1, 4)
    combined$species2 <-substr(combined[,1], 10, 13)
    category <- c() # add categories for IBD plotting (within and between)
    for(k in 1:nrow(combined)) {
      ifelse(combined[k,"species1"] == combined[k,"species2"],
             category <- c(category, paste0("within_", combined[k,"species1"])),
             category <- c(category, paste0("between_", combined[k,"species1"], "_", combined[k,"species2"]))) }
    grouping <- c() # add categories for IBD plotting (within and between)
    for(l in 1:nrow(combined)) {
      ifelse(combined[l,"species1"] == combined[l,"species2"],
             grouping <- c(grouping, paste0("within")),
             grouping <- c(grouping, paste0("between"))) }
    combined$category <- as.factor(category)
    combined$grouping <- as.factor(grouping)
    df_list[[i]][[j]] <-combined
    print(head(df_list[[i]][[j]]))
}}
save(df_list,file=paste0(sp1,"_",sp2,"_dflist.rda"))

#### plot all IBDs ####
print("plot all IBDs")
pdf(paste0("IBD_",sp1,"_",sp2,"_quickplots.pdf"),12,12)
par(mfrow=c(2,2))
par(mar=c(4.1,4.1,1.1,1.1))
  varvec<-vector()
  dvec<-vector()
  nvec<-vector()
  catvec<-vector()
  mantelps<-vector()
  mantelrs<-vector()
for(i in 1:length(data_list)){
  for(j in 1:length(distm)){
    morphod<-mat_list[[i]][[j]] ; 
    geod<-geo_list[[i]] ; 
    combined<-df_list[[i]][[j]]
    plot(log(geod+1),morphod,type="n",xlab="Geographic distance [log(km)]",
         ylab=paste0(data_types[i],"_",distm[j],"_morphodist"))
    abline(lm(morphod~log(geod+1)),lty=1,lwd=2,col=1)   
    ibd <- mantel.randtest(morphod,geod,nrepet = 999)
    mtext(substitute(italic(R)^2==rr,list(rr=round((ibd$obs)^2,3))),side=3,line=-1.1,adj=0,font = 1, cex=0.7)
    mtext(substitute(italic(p)<=pp,list(pp=round(ibd$pvalue,4))),side=3,line=-1.1,adj=1,font = 1, cex=0.7)
      # store stats
      varvec<-append(varvec,data_types[i])
      dvec<-append(dvec,distm[j])
      nvec<-append(nvec,length(morphod))
      catvec<-append(catvec,"all")
      mantelrs<-append(mantelrs,ibd$obs)
      mantelps<-append(mantelps,ibd$pvalue)
    for(k in 1:length(levels(combined$category))){
      tempdata<-combined[combined$category==levels(combined$category)[k],
                        c(which(colnames(combined)=="indPairs"),
                          which(colnames(combined)=="geoD"),
                          which(colnames(combined)=="morphoD"))]
      tempdata[,4]<-substr(tempdata[,1], 1, 8) ; tempdata[,5]<-substr(tempdata[,1], 10, 17)
        morfodist <- as.dist(t(daply(tempdata, .(V4, V5), function(x) x$morphoD)))
      print(c(i,j,k))
      print(dim(as.matrix(morfodist)))
        geodist <- as.dist(t(daply(tempdata, .(V4, V5), function(x) x$geoD)))
      print(dim(as.matrix(geodist)))
      points(log(geodist+1),morfodist,col=alpha(colorsFor2[k],0.4))
      if(length(which(geodist>0)>0)){
        print("plotting relationship")
      abline(lm(morfodist~log(geodist+1)), lty=1,lwd=2,col=colorsFor2[k])  
        ibd <- mantel.randtest(morfodist,log(geodist+1),nrepet = 999)
      legend("bottomright",legend=levels(combined$category),pch=1, col=colorsFor2,cex=0.7)  
      mtext(substitute(italic(R)^2==rr,list(rr=round((ibd$obs)^2,3))),side=3,line=-(k+1.1),adj=0,font = 1, cex=0.7,col=colorsFor2[k])
      mtext(substitute(italic(p)<=pp,list(pp=round(ibd$pvalue,4))),side=3,line=-(k+1.1),adj=1,font = 1, cex=0.7,col=colorsFor2[k])}
        varvec<-append(varvec,data_types[i])
        dvec<-append(dvec,distm[j])
        nvec<-append(nvec,length(tempdata[,1]))
        catvec<-append(catvec,levels(combined$category)[k])
        mantelrs<-append(mantelrs,ibd$obs)
        mantelps<-append(mantelps,ibd$pvalue)
    }}}
dev.off()

stat_sums<-as.data.frame(cbind(varvec,dvec,nvec,catvec,mantelps,mantelrs))
stat_sums$nvec<-as.numeric(nvec)
stat_sums$mantelrs<-as.numeric(mantelrs)
stat_sums$mantelps<-as.numeric(mantelps)

save(stat_sums,file=paste0(sp1,"_",sp2,"_stat_sums.rda"))

#### plot a few stats .... ####
print("plot a few stats")
head(stat_sums)
summary(stat_sums)
plot(stat_sums$nvec,stat_sums$mantelrs,type="p")
points(stat_sums[stat_sums$catvec==paste0("between_",sp1,"_",sp2),c("nvec","mantelrs")],col="red")
# abline(lm(stat_sums$mantelrs~stat_sums$nvec),lty=1,lwd=2,col=1) 


# which are the points with the best R stats ???
stat_sums[stat_sums$mantelrs>0.28,]


#### resampling procedure #### 
print("resampling procedure")
# lets work with the largest number of variables, because it allows for more resampling opportunities.
# let's also use the best fitting distance metric

i<-which(data_types=="12var_30NA")
j<-which(distm==stat_sums$dvec[which(stat_sums$mantelrs==max(stat_sums[stat_sums$varvec==data_types[i],"mantelrs"]))])

# get a look at all resampling possibilities and select the 
# aimed nb of resampling:
nrs<-220
samptab<-data.frame(matrix(ncol=2,nrow=last_var[i]-1))
for(n in 2:(last_var[i]-1)){
  sampcomb<-combn(colnames(data_list[[i]][,2:last_var[i]]),n)
  print("number of variables, nb of resampling options")
  print(dim(sampcomb))
  samptab[n,]<-dim(sampcomb)
}
head(sampcomb)
samptab
# let's use the max number of resamplings variables allowing for ~nrs resampling (defined above)
rsize<-max(which(abs(samptab[,2] - nrs) == min(abs(samptab[,2] - nrs), na.rm = TRUE)))
print(rsize)
resampvar<-"variables"

if(min(abs(samptab[,2] - nrs), na.rm = TRUE)>170){
  print("not enough variables shared across candidates for a proper resampling > 50")
  print("conducting resampling of individuals instead of variables")
  # get remaining number of samples in each sp candidates
  min_ind<-min(dim(data_list[[i]][data_list[[i]]$population==sp1,])[1],dim(data_list[[i]][data_list[[i]]$population==sp2,])[1])
  spmin<-which(c(dim(data_list[[i]][data_list[[i]]$population==sp1,])[1],dim(data_list[[i]][data_list[[i]]$population==sp2,])[1])==min_ind)
  samptab<-data.frame(matrix(ncol=2,nrow=min_ind-1))
  for(n in 2:(min_ind-1)){
    sampcomb<-combn(rownames(data_list[[i]][data_list[[i]]$population==c(sp1,sp2)[spmin],]),n)
    print("number of individuals, nb of resampling options")
    print(dim(sampcomb))
    samptab[n,]<-dim(sampcomb)
  }
  samptab #head(sampcomb)
  # use the max number of resamplings variables allowing for ~nrs resampling (defined above)
  rsize<-max(which(abs(samptab[,2] - nrs) == min(abs(samptab[,2] - nrs), na.rm = TRUE)))
  print(rsize)
  resampvar<-"individuals"
}

#### start resampling variables ####
if(resampvar=="variables"){
  print("resampling variables")
sampcomb<-combn(colnames(data_list[[i]][,3:last_var[i]]),rsize)
resamp<-list() # this takes ~2 sec / iteration

for(samp in 1:length(sampcomb[1,])){
  #cols<-sample(x=seq(2:13),size=rsize,replace = FALSE)
  cols<-sampcomb[,samp]
  redat<-data_list[[i]][,c("population",cols,"latitude","longitude")]
  redat$consecutive_numbers <- 100 + 1:nrow(redat) 
  redat$spp_number <- paste(redat$population, redat$consecutive_numbers, sep="_") 
  remat<-dist(redat[,2:(rsize+1)],method=distm[j])
  remat<-dist_setNames(remat, redat$spp_number)
  regeo<-as.dist(spDists(as.matrix(redat[,c("longitude","latitude")]),longlat = T))

  indPairs<-c()
  morphoD <- c()
  geoD <- c()
  morphomat<-as.matrix(remat)
  geomat<-as.matrix(regeo)
  for(row in 1:nrow(morphomat)){
    for(col in 1:ncol(morphomat)){
      if (row < col) {
        indPairs <- c(indPairs,paste0(rownames(morphomat)[row],"_",
                                      colnames(morphomat)[col]))
        morphoD <- c(morphoD, morphomat[row, col])
        geoD <- c(geoD, geomat[row, col])
      }}}
  combined <- data.frame(indPairs, morphoD, geoD)
  combined$species1 <-substr(combined[,1], 1, 4)
  combined$species2 <-substr(combined[,1], 10, 13)
  category <- c() # add categories for IBD plotting (within and between)
  for(k in 1:nrow(combined)) {
    ifelse(combined[k,"species1"] == combined[k,"species2"],
           category <- c(category, paste0("within_", combined[k,"species1"])),
           category <- c(category, paste0("between_", combined[k,"species1"], "_", combined[k,"species2"]))) }
  grouping <- c() # add categories for IBD plotting (within and between)
  for(l in 1:nrow(combined)) {
    ifelse(combined[l,"species1"] == combined[l,"species2"],
           grouping <- c(grouping, paste0("within")),
           grouping <- c(grouping, paste0("between"))) }
  combined$category <- as.factor(category)
  combined$grouping <- as.factor(grouping)
  resamp[[samp]]<-list()
  resamp[[samp]][[1]] <-redat
  resamp[[samp]][[2]] <-remat
  resamp[[samp]][[3]] <-regeo
  resamp[[samp]][[4]] <-combined
  print(samp)
  print(head(combined))
}}

#### start resampling individuals instead ####
if(resampvar=="individuals"){
  print("resampling individuals")
  data_list[[i]]<-as.data.frame(data_list[[i]])
  if(spmin==1){
    sampcomb1<-combn(rownames(data_list[[i]][data_list[[i]]$population==sp1,]),rsize)
    sampcomb2<-sampcomb1
    for(nbrsize in 1:length(sampcomb2[1,])){
      sampcomb2[,nbrsize]<-sample(rownames(data_list[[i]][data_list[[i]]$population==sp2,]),rsize)}
    sampcomb<-rbind(sampcomb1,sampcomb2)}
  if(spmin==2){
    sampcomb1<-combn(rownames(data_list[[i]][data_list[[i]]$population==sp2,]),rsize)
    sampcomb2<-sampcomb1
    for(nbrsize in 1:length(sampcomb2[1,])){
      sampcomb2[,nbrsize]<-sample(rownames(data_list[[i]][data_list[[i]]$population==sp1,]),rsize)}
    sampcomb<-rbind(sampcomb1,sampcomb2)}
  #dim(sampcomb1)
  #dim(sampcomb2)
  dim(sampcomb)
  #sampcomb1[1:rsize,31:35] ; sampcomb2[1:rsize,31:35] 
  
  resamp<-list() # this takes ~2 sec / iteration
  
  for(samp in 1:length(sampcomb[1,])){
    #cols<-sample(x=seq(2:13),size=rsize,replace = FALSE)
    cols<-sampcomb[,samp]
    redat<-data_list[[i]][cols,]
    redat$consecutive_numbers <- 100 + 1:nrow(redat) 
    redat$spp_number <- paste(redat$population, redat$consecutive_numbers, sep="_") 
    remat<-dist(redat[,2:last_var[i]],method=distm[j])
    remat<-dist_setNames(remat, redat$spp_number)
    regeo<-as.dist(spDists(as.matrix(redat[,c("longitude","latitude")]),longlat = T))
    
    indPairs<-c()
    morphoD <- c()
    geoD <- c()
    morphomat<-as.matrix(remat)
    geomat<-as.matrix(regeo)
    for(row in 1:nrow(morphomat)){
      for(col in 1:ncol(morphomat)){
        if (row < col) {
          indPairs <- c(indPairs,paste0(rownames(morphomat)[row],"_",
                                        colnames(morphomat)[col]))
          morphoD <- c(morphoD, morphomat[row, col])
          geoD <- c(geoD, geomat[row, col])
        }}}
    combined <- data.frame(indPairs, morphoD, geoD)
    combined$species1 <-substr(combined[,1], 1, 4)
    combined$species2 <-substr(combined[,1], 10, 13)
    category <- c() # add categories for IBD plotting (within and between)
    for(k in 1:nrow(combined)) {
      ifelse(combined[k,"species1"] == combined[k,"species2"],
             category <- c(category, paste0("within_", combined[k,"species1"])),
             category <- c(category, paste0("between_", combined[k,"species1"], "_", combined[k,"species2"]))) }
    grouping <- c() # add categories for IBD plotting (within and between)
    for(l in 1:nrow(combined)) {
      ifelse(combined[l,"species1"] == combined[l,"species2"],
             grouping <- c(grouping, paste0("within")),
             grouping <- c(grouping, paste0("between"))) }
    combined$category <- as.factor(category)
    combined$grouping <- as.factor(grouping)
    resamp[[samp]]<-list()
    resamp[[samp]][[1]] <-redat
    resamp[[samp]][[2]] <-remat
    resamp[[samp]][[3]] <-regeo
    resamp[[samp]][[4]] <-combined
    print(samp)
    print(head(combined))
  }}

save(resamp,file=paste0(sp1,"_",sp2,"_resamp.rda"))

#### extract IBD stats from resampled data ####
varvecs<-vector()
dvecs<-vector()
nvecs<-vector()
catvecs<-vector()
mantelpss<-vector()
mantelrss<-vector()
for(samp in 1:length(resamp)){
  #for(j in 1:length(distm)){
    morphod<-resamp[[samp]][[2]] ; 
    geod<-resamp[[samp]][[3]] ; 
    combined<-resamp[[samp]][[4]]
    #plot(log(geod+1),morphod,type="n",xlab="Geographic distance [log(km)]",
    #     ylab=paste0(data_types[i],"_",distm[j],"_morphodist"))
    #abline(lm(morphod~log(geod+1)),lty=1,lwd=2,col=1)   
    ibd <- mantel.randtest(morphod,geod,nrepet = 999)
    #mtext(substitute(italic(R)^2==rr,list(rr=round((ibd$obs)^2,3))),side=3,line=-1.1,adj=0,font = 1, cex=0.7)
    #mtext(substitute(italic(p)<=pp,list(pp=round(ibd$pvalue,4))),side=3,line=-1.1,adj=1,font = 1, cex=0.7)
    # store stats
    varvecs<-append(varvecs,data_types[i])
    dvecs<-append(dvecs,distm[j])
    nvecs<-append(nvecs,length(morphod))
    catvecs<-append(catvecs,"all")
    mantelrss<-append(mantelrss,ibd$obs)
    mantelpss<-append(mantelpss,ibd$pvalue)
    for(k in 1:length(levels(combined$category))){
      tempdata<-combined[combined$category==levels(combined$category)[k],
                         c(which(colnames(combined)=="indPairs"),
                           which(colnames(combined)=="geoD"),
                           which(colnames(combined)=="morphoD"))]
      tempdata[,4]<-substr(tempdata[,1], 1, 8) ; tempdata[,5]<-substr(tempdata[,1], 10, 17)
      morfodist <- as.dist(t(daply(tempdata, .(V4, V5), function(x) x$morphoD)))
      print(c(samp,i,j,k))
      print(dim(as.matrix(morfodist)))
      geodist <- as.dist(t(daply(tempdata, .(V4, V5), function(x) x$geoD)))
      print(dim(as.matrix(geodist)))
      #points(log(geodist+1),morfodist,col=alpha(colorsFor2[k],0.4))
      #abline(lm(morfodist~log(geodist+1)), lty=1,lwd=2,col=colorsFor2[k])  
      ibd <- mantel.randtest(morfodist,log(geodist+1),nrepet = 999)
      #legend("bottomright",legend=levels(combined$category),pch=1, col=colorsFor2,cex=0.7)  
      #mtext(substitute(italic(R)^2==rr,list(rr=round((ibd$obs)^2,3))),side=3,line=-(k+1.1),adj=0,font = 1, cex=0.7,col=colorsFor2[k])
      #mtext(substitute(italic(p)<=pp,list(pp=round(ibd$pvalue,4))),side=3,line=-(k+1.1),adj=1,font = 1, cex=0.7,col=colorsFor2[k])
      varvecs<-append(varvecs,data_types[i])
      dvecs<-append(dvecs,distm[j])
      nvecs<-append(nvecs,length(tempdata[,1]))
      catvecs<-append(catvecs,levels(combined$category)[k])
      mantelrss<-append(mantelrss,ibd$obs)
      mantelpss<-append(mantelpss,ibd$pvalue)
    }}#}


stat_samp<-as.data.frame(cbind(varvecs,dvecs,nvecs,catvecs,mantelpss,mantelrss))
stat_samp$nvecs<-as.numeric(nvecs)
stat_samp$mantelrss<-as.numeric(mantelrss)
stat_samp$mantelpss<-as.numeric(mantelpss)

save(stat_samp,file=paste0(sp1,"_",sp2,"_stat_samp.rda"))

#### plot a few stats ####
head(stat_samp)
summary(stat_samp)
plot(stat_samp$nvecs,stat_samp$mantelrss,type="p")
points(stat_samp[stat_samp$catvecs=="between_arno_Msp1",c("nvecs","mantelrss")],col="red")
#abline(lm(stat_samp$mantelrs~stat_samp$nvecs),lty=1,lwd=2,col=1) 

ggplot(stat_samp,aes(x=mantelrss,fill=catvecs))+
  geom_density(alpha=.3)

hist(stat_samp$mantelrss,nclass=40)

#### estimate lognrmse ####
source("./createLOGNrmse_function.R")
s_within<-vector()
s_between<-vector()
s_nrmse<-vector()
rm(mydf)
for(samp in 1:length(resamp)){
  mydf<-resamp[[samp]][[4]]
  print(paste0("processing sample : ",samp)) #head(mydf)
  mydf$combo<-paste0(mydf$species1,"-",mydf$species2)
  # create log(geo) column with epsilon correction # 
  c_epsilon = as.numeric(quantile(mydf$geoD,na.rm = T)[2])
  if(c_epsilon==0){
    c_epsilon = 0.01
  }else{}
  mydf$log_geo_dist = log(mydf$geoD + c_epsilon)
  mydf$value<-mydf$morphoD
  mylist_pairs = unique(mydf$combo)
  lgn<-LOGnrme_regression(mydf,mylist_pairs)
  s_within<-append(s_within,as.character(lgn$within))
  s_between<-append(s_between,as.character(lgn$between))
  s_nrmse<-append(s_nrmse,lgn[,3])
}

lrmse_samp<-as.data.frame(cbind(s_within,s_between,s_nrmse))
lrmse_samp$s_nrmse<-as.numeric(s_nrmse)
summary(lrmse_samp)

save(lrmse_samp,file=paste0(sp1,"_",sp2,"_lrmse_samp.rda"))




#### plot lognmrse ####
pdf(paste0("IBD_",sp1,"_",sp2,"_lognrmseplots.pdf"),12,12)
myplot<-ggplot(lrmse_samp,aes(x=s_nrmse,fill=s_within))+
  geom_density(alpha=.3)+
  xlim(0,1)
print(myplot)
dev.off()

#### end of the species pair loop ####
}}}

#### Plot all lognrmse together #### 
lrmse_samp_all<-data.frame(matrix(ncol=3,nrow=0))
lrmse_files<-list.files(path=".",pattern = "_lrmse_samp.rda")
for(pair in 1:length(lrmse_files)){
  load(lrmse_files[pair],verbose=T)
  lrmse_samp_all<-rbind(lrmse_samp_all,lrmse_samp)}
#lrmse_samp_all$s_within<-as.factor(lrmse_samp_all$s_within)
lrmse_samp_all[,1]<-substr(lrmse_samp_all[,1], 1, 4)
lrmse_samp_all[,1]<-as.factor(lrmse_samp_all[,1])
summary(lrmse_samp_all)
allplot<-ggplot(lrmse_samp_all,aes(x=s_nrmse,fill=s_within))+
  geom_density(alpha=.3)+
  #geom_vline(xintercept = quantile(lrmse_samp_all[lrmse_samp_all[,1]=="mitt",3],probs = 0.95),lwd=0.7,colour="blue")+
  #geom_vline(xintercept = quantile(lrmse_samp_all[lrmse_samp_all[,1]=="lehi",3],probs = 0.95),lwd=0.7,colour="green")+
  #geom_vline(xintercept = quantile(lrmse_samp_all[lrmse_samp_all[,1]=="tavN" & lrmse_samp_all[,2]=="tavN-tavC" ,3],probs = 0.95),lwd=0.7,colour="pink")+
  geom_vline(xintercept = quantile(lrmse_samp_all[lrmse_samp_all[,2]=="tavN-tavC" ,3],probs = 0.95),lwd=0.7,colour="red")+
  facet_grid(rows=vars(s_between),scales = "free_y")+
  xlim(0.05,0.7) + theme(legend.position = "none") +
  xlab("Log[nRMSE]") + ylab("")

allplot

pdf(paste0("IBD_all_lognrmseplots.pdf"),9,12)
print(allplot)
dev.off()



#### Plot all pvals and R² together #### 
stat_samp_all<-data.frame(matrix(ncol=6,nrow=0))
stat_files<-list.files(path=".",pattern = "_stat_samp.rda")
for(pair in 1:length(stat_files)){
  load(stat_files[pair],verbose=T)
  stat_samp$pair<-rep(paste0(substr(stat_files[pair], 1, 4),"-",substr(stat_files[pair], 6, 9)),times=length(stat_samp[,1]))
  stat_samp<-stat_samp[stat_samp$catvecs!=paste0("between","_",substr(stat_files[pair], 1, 4),"_",substr(stat_files[pair], 6, 9)),]
  stat_samp_all<-rbind(stat_samp_all,stat_samp)}
stat_samp_all$pair<-as.factor(stat_samp_all$pair)
stat_samp_all <-na.omit(stat_samp_all)
stat_samp_all <- stat_samp_all[stat_samp_all$pair!="NA-NA",]
stat_samp_all<-droplevels(stat_samp_all[stat_samp_all$catvecs=="all",])
# change twisted level orders if there are some
levels(stat_samp_all$pair)[levels(stat_samp_all$pair)=="gerp-joll"] <- "joll-gerp"

summary(stat_samp_all)
allpplot<-ggplot(stat_samp_all,aes(x=mantelpss,fill=catvecs))+
  geom_density(alpha=.3)+
  geom_vline(xintercept = 0.05,lwd=0.7,colour="red")+
  facet_grid(rows=vars(pair),scales = "free_y")+
  xlim(0,1) + scale_x_log10() +
  theme(legend.position = "none") +
  xlab("Mantel p-value") + ylab("Density")
allpplot

pdf(paste0("IBD_all_pvals_plots.pdf"),9,12)
print(allpplot)
dev.off()

allrplot<-ggplot(stat_samp_all,aes(x=mantelrss,fill=catvecs))+
  geom_density(alpha=.3)+
  facet_grid(rows=vars(pair),scales = "free_y")+
  xlim(0,1) + #scale_x_log10() +
  theme(legend.position = "none") +
  #geom_vline(xintercept = quantile(stat_samp_all[stat_samp_all$pair=="mitt-lehi",6],probs = 0.05),lwd=0.7,colour="green")+
  geom_vline(xintercept = quantile(stat_samp_all[stat_samp_all$pair=="tavN-tavC" ,6],probs = 0.05),lwd=0.7,colour="red")+
  xlab("Mantel R") + ylab("")
allrplot
pdf(paste0("IBD_all_rvals_plots.pdf"),9,12)
print(allrplot)
dev.off()

# combine plots
figure <- ggarrange(ggarrange(allpplot,allrplot,labels=c("A","B"),ncol=2,nrow = 1),allplot,labels=c("","C"),ncol=2,nrow = 1)

pdf(paste0("IBD_all_3stat_plots.pdf"),12,12)
print(figure)
dev.off()


#### Plot lognrmse and p for each pair and triplet separately #### 
pairs<-list.dirs(path="D:/jordi/mada/consortium_rad_microcebus/js_phylo/ngsAdmix",full.names = FALSE,recursive = FALSE)
pcalist<-list.files(path="D:/jordi/mada/consortium_rad_microcebus/js_phylo/IBD/IBD_morpho/pca_morpho/",pattern = ".rda")
pcalist<-gsub(pattern = ".rda","",pcalist)
pcalist<-pcalist[c(1:5,7,6)]
mouse_colors<-read.csv(file = "infoColors_all_mouse_lemurs.csv")
pairs
#pairs2<-pairs
pairs2<-gsub("_","",pairs)
pairs2x<-pairs2[pairs2!="borsim" & pairs2!="ganmanmur"]
lrmse_files<-list.files(path=".",pattern = "_lrmse_samp.rda")
niv<-levels(lrmse_samp_all$s_between)
niv2<-gsub("_lrmse_samp.rda","",lrmse_files)
groups<-list()
groups[[1]]<-niv2[1]; groups[[2]]<-niv2[c(2,11,12)]; groups[[3]]<-niv2[c(3,13,14)]; groups[[4]]<-niv2[4]; 
groups[[5]]<-niv2[5]; groups[[6]]<-niv2[10]; groups[[7]]<-niv2[c(6,7,9)];

for(tri in 1:length(groups)){
  lrmse_samp_g<-data.frame(matrix(ncol=3,nrow=0))
  stat_samp_g<-data.frame(matrix(ncol=6,nrow=0))
  if(length(groups[[tri]])==1){
    lrmse_files<-list.files(path=".",pattern = paste0(groups[[tri]][1],"_lrmse_samp.rda"))}
  if(length(groups[[tri]])==2){
    lrmse_files<-list.files(path=".",pattern = c(paste0(groups[[tri]][1],"_lrmse_samp.rda"),
                                                 paste0(groups[[tri]][2],"_lrmse_samp.rda")))}
  if(length(groups[[tri]])==3){
    lrmse_files<-list.files(path=".",pattern = c(paste0(groups[[tri]][1],"_lrmse_samp.rda"),
                                                 paste0(groups[[tri]][2],"_lrmse_samp.rda"),
                                                 paste0(groups[[tri]][3],"_lrmse_samp.rda")))}
  for(pair in 1:length(groups[[tri]])){
  load(paste0(groups[[tri]][pair],"_lrmse_samp.rda"),verbose=T)
  lrmse_samp[,1]<-substr(lrmse_samp[,1], 1, 4)
  lrmse_samp_g<-rbind(lrmse_samp_g,lrmse_samp)}
  lrmse_samp_g[,1]<-substr(lrmse_samp_g[,1], 1, 4)
  lrmse_samp_g[,1]<-as.factor(lrmse_samp_g[,1])
  summary(lrmse_samp_g)
  
mycolors<-as.character(droplevels(mouse_colors$code[mouse_colors$sp_code2 %in% levels(lrmse_samp_g$s_within)]))
  
  rmseplot<-ggplot(lrmse_samp_g,aes(x=s_nrmse,fill=s_within,color=s_within))+
  geom_density(aes(fill = s_within),alpha=.5)+
  scale_color_manual(name='',values = mycolors)+
  scale_fill_manual(name='',labels = names(mycolors), values=mycolors)+
  geom_vline(xintercept = quantile(lrmse_samp_all[lrmse_samp_all[,2]=="tavN-tavC" ,3],probs = 0.95),lwd=0.7,colour="red")+
  facet_grid(rows=vars(s_between),scales = "free_y")+
  xlab("Log[nRMSE]") + ylab("")+ theme_bw() + theme(legend.position = "none",strip.text.x = element_text(size=2)) + 
    #theme(legend.position = c(0.1, 0.9),legend.title= element_blank())+
  scale_x_log10(limits = c(0.09, 0.9), expand = c(0, 0)) #xlim(0.05,0.7) + 
  rmseplot <- rmseplot + theme(strip.text = element_text(size = 6))
  rmseplot
  
## pvals
  for(pair in 1:length(groups[[tri]])){
  load(paste0(groups[[tri]][pair],"_stat_samp.rda"),verbose=T)
 
  stat_samp$pair<-rep(paste0(substr(groups[[tri]][pair], 1, 4),"-",substr(groups[[tri]][pair], 6, 9)),times=length(stat_samp[,1]))
  stat_samp<-stat_samp[stat_samp$catvecs!=paste0("between","_",substr(groups[[tri]][pair], 1, 4),"_",substr(groups[[tri]][pair], 6, 9)),]
  stat_samp_g<-rbind(stat_samp_g,stat_samp)}
  stat_samp_g$pair<-as.factor(stat_samp_g$pair)
  stat_samp_g <-na.omit(stat_samp_g)
  stat_samp_g <- stat_samp_g[stat_samp_g$pair!="NA-NA",]
  stat_samp_g<-droplevels(stat_samp_g[stat_samp_g$catvecs=="all",])
  summary(stat_samp_g)
  #stat_samp$pair<-rep(paste0(substr(stat_files[pair], 1, 4),"-",substr(stat_files[pair], 6, 9)),times=length(stat_samp[,1]))
  #stat_samp<-stat_samp[stat_samp$catvecs!=paste0("between","_",substr(stat_files[pair], 1, 4),"_",substr(stat_files[pair], 6, 9)),]
  #stat_samp<-droplevels(stat_samp[stat_samp$catvecs=="all",])
pplot<-ggplot(stat_samp_g,aes(x=mantelpss))+
  geom_density(alpha=.3,color="darkblue", fill="lightblue")+
  #geom_density_ridges(aes(fill = catvecs),alpha =.7,scale =.9,size=.3,show.legend = F)+
  geom_vline(xintercept = 0.05,lwd=0.7,colour="red")+
  facet_grid(rows=vars(pair),scales = "free_y")+
  theme(legend.position = "none") + scale_x_log10(limits = c(0.001, 1), expand = c(0, 0))+
  xlab("p-value") + ylab("density") + theme_bw()
pplot <- pplot + theme(strip.text = element_text(size = 6))
pplot

# load pca plot
load(paste0("./pca_morpho/",pcalist[tri],".rda"),verbose=T)


## combine figures
print("combining figure")
pairfigure <- ggarrange(pplot,rmseplot,ncol=2,nrow = 1)
pairpcafigure <- ggarrange(pairfigure,pcagplot,heights=c(1,2),ncol=1,nrow = 2)
#png(paste0("IBD_lognrmseplots_",pairs2x[tri],".3.png"), width = 6.35, height = 2.6, units = 'in', res = 300)
png(paste0("IBD_lognrmseplots_",pairs2x[tri],".4.png"), width = 6.35, height = 8.67, units = 'in', res = 300)

print(pairpcafigure)
dev.off()
}


#### Build result table ####
taxo_tab<-as.data.frame(matrix(nrow=length(levels(lrmse_samp_all$s_between)),ncol=11))
taxo_tab$V1<-levels(lrmse_samp_all$s_between)
taxo_tab
for(i in 1:length(taxo_tab[,1])){
  # pvals
  taxo_tab[i,2]<-quantile(stat_samp_all[stat_samp_all$pair==taxo_tab[i,1],5],probs = 0.95)
  # rvals 05 - 95
  taxo_tab[i,3:4]<-quantile(stat_samp_all[stat_samp_all$pair==taxo_tab[i,1],6],probs = c(0.05,0.95))
  # log-nrmse 05 - 95 
  taxo_tab[i,5:7]<-quantile(lrmse_samp_all[lrmse_samp_all$s_between==taxo_tab[i,1] & 
                            lrmse_samp_all$s_within==strsplit(taxo_tab[i,1],"-")[[1]][1],3],probs = c(0.05,0.5,0.95))
  taxo_tab[i,8:10]<-quantile(lrmse_samp_all[lrmse_samp_all$s_between==taxo_tab[i,1] & 
                                             lrmse_samp_all$s_within==strsplit(taxo_tab[i,1],"-")[[1]][2],3],probs = c(0.05,0.5,0.95))
  if(taxo_tab[i,2]>0.05){
    taxo_tab[i,11]<-"non-significant IBD"
  } else if (taxo_tab[i,4]<quantile(stat_samp_all[stat_samp_all$pair=="tavN-tavC",6],probs = 0.05)){
      taxo_tab[i,11]<-"low IBD - R<q05-tavN-tavC"
  }  else if (taxo_tab[i,7]<=quantile(lrmse_samp_all[lrmse_samp_all$s_between=="tavN-tavC" & 
                                                    lrmse_samp_all$s_within=="tavN",3],probs = 0.95) | 
             taxo_tab[i,10]<=quantile(lrmse_samp_all[lrmse_samp_all$s_between=="tavN-tavC" & 
                                                   lrmse_samp_all$s_within=="tavN",3],probs = 0.95)){
      taxo_tab[i,11]<-"clear morpho IBD"
  }  else if (taxo_tab[i,5]>quantile(lrmse_samp_all[lrmse_samp_all$s_between=="tavN-tavC" & 
                                                    lrmse_samp_all$s_within=="tavC",3],probs = 0.95) & 
              taxo_tab[i,8]>quantile(lrmse_samp_all[lrmse_samp_all$s_between=="tavN-tavC" & 
                                                    lrmse_samp_all$s_within=="tavC",3],probs = 0.95)){
    taxo_tab[i,11]<-"low-fit log-nrmse"
  }  else if (taxo_tab[i,6]<quantile(lrmse_samp_all[lrmse_samp_all$s_between=="tavN-tavC" & 
                                                    lrmse_samp_all$s_within=="tavC",3],probs = 0.95) | 
              taxo_tab[i,9]<quantile(lrmse_samp_all[lrmse_samp_all$s_between=="tavN-tavC" & 
                                                    lrmse_samp_all$s_within=="tavC",3],probs = 0.95)){
    taxo_tab[i,11]<-"high-fit log-nrmse [IBD]"
  } else {taxo_tab[i,11]<-"grey-zone [IBD]"}
}
colnames(taxo_tab)<-c("pair","pval-q95","R-q05","R-q95","p1lnrmse-q05","p1lnrmse-q50","p1lnrmse-q95","p2lnrmse-q05","p2lnrmse-q50","p2lnrmse-q95","results")
taxo_tab
write.table(taxo_tab,file="taxo_table.tsv",quote = F,sep="\t",row.names = F,col.names = T)  
  
head(lrmse_samp_all)
head(stat_samp_all)

#### Plot all IBD together #### 
ibd_all_dat<-data.frame(matrix(ncol=3,nrow=length(lrmse_files)))
ibd_all_dat[,1]<-substr(lrmse_files, 1, 4)
ibd_all_dat[,2]<-substr(lrmse_files, 6, 9)
head(ibd_all_dat)



