### Script by Dominik Schuessler

# Calculation of species-specific ENMs with individual fine tuning of parameters.
# This is followed by the calculation of niche overlap.
# Note: ENMtools and accordingly this script heavily depends on the "raster" package;
# it can be installed with: 
# devtools::install_github("danlwarren/ENMTools", ref = "raster")

#### +++ Main script +++ ####
#### Packages ####
library(sp)
library(raster)     # handling of GIS data
library(terra)      # handling of GIS data
library(dismo)
library(rJava)
library(ENMTools)   # ENM calculation
library(ENMeval)    # parameter selection, fine tuning of ENMs
library(RStoolbox)
library(dplyr)      # re-arrange data
library(readxl)     # read Excel-table
library(ecospat)    # to calculate CBI
library(xlsx)       # write Excel-table





#### Data preparation ####
setwd("C://Users//domin//OneDrive//Dokumente//011-Revised_taxonomy_of_Microcebus//ENMs")
data <- read_excel("02_Occurrence_data_Microcebus_20230612.xlsx", sheet = "data", col_names = TRUE)
set.seed(21)
options(max.print=100000)

length(unique(data$species)) # number of taxa in dataset




### load bioclimatic data (CHELSA 2.1 database)
bio03 <- rast("C://Users//domin//Documents//999_Bioclim_Mada//Chelsa_2.1//CH_bio03.tif")
bio04 <- rast("C://Users//domin//Documents//999_Bioclim_Mada//Chelsa_2.1//CH_bio04.tif")
bio05 <- rast("C://Users//domin//Documents//999_Bioclim_Mada//Chelsa_2.1//CH_bio05.tif")
bio06 <- rast("C://Users//domin//Documents//999_Bioclim_Mada//Chelsa_2.1//CH_bio06.tif")
bio12 <- rast("C://Users//domin//Documents//999_Bioclim_Mada//Chelsa_2.1//CH_bio12.tif")
bio15 <- rast("C://Users//domin//Documents//999_Bioclim_Mada//Chelsa_2.1//CH_bio15.tif")
bio16 <- rast("C://Users//domin//Documents//999_Bioclim_Mada//Chelsa_2.1//CH_bio16.tif")
bio17 <- rast("C://Users//domin//Documents//999_Bioclim_Mada//Chelsa_2.1//CH_bio17.tif")

# stack bioclimatic layers and perform PCA
climate_variables_raw <- c(bio03, bio04, bio05, bio06, bio12, bio15, bio16, bio17) # stack all layers into one dataset
pca1 <- rasterPCA(climate_variables_raw, spca = TRUE) # calculate PCA of bioclimatic variables
summary(pca1$model)  # get model summary (i.e., proportion of variance of each PC)
loadings(pca1$model) # get loading on PCs

# extract first 3 PCs, which explain together >93% of the variation in the data set
PC1 <- pca1$map$PC1
PC2 <- pca1$map$PC2
PC3 <- pca1$map$PC3

climate_variables <- c(PC1, PC2, PC3) # stack first 3 PCs
plot(climate_variables) # visualize to check data validity




### preparation of species data sets 
str(data) # data overview
data$long <- as.numeric(data$long) # transformation to numeric data structure
data$lat <- as.numeric(data$lat)   # transformation to numeric data structure
reduced_df <- dplyr::select(data, species, long, lat) # downsize data frame
reduced_df <- na.omit(reduced_df) # drop all rows with missing data

# reduce to taxon subsets 
df_rave <- filter(reduced_df, species == "Mravelobensis")
df_bong <- filter(reduced_df, species == "Mbongolavensis")
df_danf <- filter(reduced_df, species == "Mdanfossi")

df_myox <- filter(reduced_df, species == "Mmyoxinus")
df_bert <- filter(reduced_df, species == "Mberthae")
df_rufu <- filter(reduced_df, species == "Mrufus")

df_lehi <- filter(reduced_df, species == "Mlehilahytsara")
df_mitt <- filter(reduced_df, species == "Mmittermeieri")

df_muriN <- filter(reduced_df, species == "MmurinusN")
df_muriC <- filter(reduced_df, species == "MmurinusC")
df_muriS <- filter(reduced_df, species == "Mganzhorni") # is all murinus like taxa in SE
df_gris <- filter(reduced_df, species == "Mgriseorufus")

df_maca <- filter(reduced_df, species == "Mmacarthurii")
df_jona <- filter(reduced_df, species == "Mjonahi")

df_arnh <- filter(reduced_df, species == "Marnholdi")
df_sp1 <- filter(reduced_df, species == "Msp1")

df_tava <- filter(reduced_df, species == "Mtavaratra")

df_tano <- filter(reduced_df, species == "Mtanosi")

df_samb <- filter(reduced_df, species == "Msambiranensis")
df_marg <- filter(reduced_df, species == "Mmargotmarshae")
df_mami <- filter(reduced_df, species == "Mmamiratra")

df_gerp <- filter(reduced_df, species == "Mgerpi")
df_joll <- filter(reduced_df, species == "Mjollyae")
df_maro <- filter(reduced_df, species == "Mmarohita")

df_simm <- filter(reduced_df, species == "Msimmonsi")
df_bora <- filter(reduced_df, species == "Mboraha")


# convert to SpatVector-object with WGS1984 geographic coordinate system (CRS)
rave <- vect(df_rave, geom = c("long", "lat"), crs = "EPSG:4326")
bong <- vect(df_bong, geom = c("long", "lat"), crs = "EPSG:4326")
danf <- vect(df_danf, geom = c("long", "lat"), crs = "EPSG:4326")

myox <- vect(df_myox, geom = c("long", "lat"), crs = "EPSG:4326")
bert <- vect(df_bert, geom = c("long", "lat"), crs = "EPSG:4326")
rufu <- vect(df_rufu, geom = c("long", "lat"), crs = "EPSG:4326")

lehi <- vect(df_lehi, geom = c("long", "lat"), crs = "EPSG:4326")
mitt <- vect(df_mitt, geom = c("long", "lat"), crs = "EPSG:4326")

muriN <- vect(df_muriN, geom = c("long", "lat"), crs = "EPSG:4326")
muriC <- vect(df_muriC, geom = c("long", "lat"), crs = "EPSG:4326")
muriS <- vect(df_muriS, geom = c("long", "lat"), crs = "EPSG:4326") # is all murinus like taxa in SE
gris <- vect(df_gris, geom = c("long", "lat"), crs = "EPSG:4326") 

maca <- vect(df_maca, geom = c("long", "lat"), crs = "EPSG:4326")
jona <- vect(df_jona, geom = c("long", "lat"), crs = "EPSG:4326")

arnh <- vect(df_arnh, geom = c("long", "lat"), crs = "EPSG:4326")
sp1 <- vect(df_sp1, geom = c("long", "lat"), crs = "EPSG:4326")

tava <- vect(df_tava, geom = c("long", "lat"), crs = "EPSG:4326")

tano <- vect(df_tano, geom = c("long", "lat"), crs = "EPSG:4326")

samb <- vect(df_samb, geom = c("long", "lat"), crs = "EPSG:4326")
marg <- vect(df_marg, geom = c("long", "lat"), crs = "EPSG:4326")
mami <- vect(df_mami, geom = c("long", "lat"), crs = "EPSG:4326")

gerp <- vect(df_gerp, geom = c("long", "lat"), crs = "EPSG:4326")
joll <- vect(df_joll, geom = c("long", "lat"), crs = "EPSG:4326")
maro <- vect(df_maro, geom = c("long", "lat"), crs = "EPSG:4326")

simm <- vect(df_simm, geom = c("long", "lat"), crs = "EPSG:4326")
bora <- vect(df_bora, geom = c("long", "lat"), crs = "EPSG:4326")



#### Setup enmtools.species objects ####
# => to be done for each taxon separately
# => sampling background should include whole Madagascar, 
#    to include the same bioclimatic conditions for all taxa. 

### M. ravelobensis
enm.rave <- enmtools.species()
enm.rave$species.name <- "Mravelobensis"
enm.rave$presence.points <- rave
enm.rave$range <- background.raster.buffer(enm.rave$presence.points, 50000, mask = climate_variables)
enm.rave$background.points <- background.points.buffer(points = enm.rave$presence.points,
                                                           radius = 1000000, n = 10000, mask = climate_variables[[1]])
### M. danfossi
enm.danf <- enmtools.species()
enm.danf$species.name <- "Mdanfossi"
enm.danf$presence.points <- danf
enm.danf$range <- background.raster.buffer(enm.danf$presence.points, 50000, mask = climate_variables)
enm.danf$background.points <- background.points.buffer(points = enm.danf$presence.points,
                                                       radius = 1000000, n = 10000, mask = climate_variables[[1]])
### M. bongolavensis
enm.bong <- enmtools.species()
enm.bong$species.name <- "Mbongolavensis"
enm.bong$presence.points <- bong
enm.bong$range <- background.raster.buffer(enm.bong$presence.points, 50000, mask = climate_variables)
enm.bong$background.points <- background.points.buffer(points = enm.bong$presence.points,
                                                       radius = 1000000, n = 10000, mask = climate_variables[[1]])

### M. myoxinus
enm.myox <- enmtools.species()
enm.myox$species.name <- "Mmyoxinus"
enm.myox$presence.points <- myox
enm.myox$range <- background.raster.buffer(enm.myox$presence.points, 50000, mask = climate_variables)
enm.myox$background.points <- background.points.buffer(points = enm.myox$presence.points,
                                                       radius = 1000000, n = 10000, mask = climate_variables[[1]])
### M. rufus
enm.rufu <- enmtools.species()
enm.rufu$species.name <- "Mrufus"
enm.rufu$presence.points <- rufu
enm.rufu$range <- background.raster.buffer(enm.rufu$presence.points, 50000, mask = climate_variables)
enm.rufu$background.points <- background.points.buffer(points = enm.rufu$presence.points,
                                                       radius = 1000000, n = 10000, mask = climate_variables[[1]])
### M. berthae
enm.bert <- enmtools.species()
enm.bert$species.name <- "Mberthae"
enm.bert$presence.points <- bert
enm.bert$range <- background.raster.buffer(enm.bert$presence.points, 50000, mask = climate_variables)
enm.bert$background.points <- background.points.buffer(points = enm.bert$presence.points,
                                                       radius = 1000000, n = 10000, mask = climate_variables[[1]])

### M. simmonsi
enm.simm <- enmtools.species()
enm.simm$species.name <- "Msimmonsi"
enm.simm$presence.points <- simm
enm.simm$range <- background.raster.buffer(enm.simm$presence.points, 50000, mask = climate_variables)
enm.simm$background.points <- background.points.buffer(points = enm.simm$presence.points,
                                                       radius = 1000000, n = 10000, mask = climate_variables[[1]])

### M. macarthurii
enm.maca <- enmtools.species()
enm.maca$species.name <- "Mmacarthurii"
enm.maca$presence.points <- maca
enm.maca$range <- background.raster.buffer(enm.maca$presence.points, 50000, mask = climate_variables)
enm.maca$background.points <- background.points.buffer(points = enm.maca$presence.points,
                                                       radius = 1000000, n = 10000, mask = climate_variables[[1]])
### M. jonahi
enm.jona <- enmtools.species()
enm.jona$species.name <- "Mjonahi"
enm.jona$presence.points <- jona
enm.jona$range <- background.raster.buffer(enm.jona$presence.points, 50000, mask = climate_variables)
enm.jona$background.points <- background.points.buffer(points = enm.jona$presence.points,
                                                       radius = 1000000, n = 10000, mask = climate_variables[[1]])

### M. tavaratra
enm.tava <- enmtools.species()
enm.tava$species.name <- "Mtavaratra"
enm.tava$presence.points <- tava
enm.tava$range <- background.raster.buffer(enm.tava$presence.points, 50000, mask = climate_variables)
enm.tava$background.points <- background.points.buffer(points = enm.tava$presence.points,
                                                       radius = 1000000, n = 10000, mask = climate_variables[[1]])

### M. sp1
enm.sp1 <- enmtools.species()
enm.sp1$species.name <- "Msp1"
enm.sp1$presence.points <- sp1
enm.sp1$range <- background.raster.buffer(enm.sp1$presence.points, 50000, mask = climate_variables)
enm.sp1$background.points <- background.points.buffer(points = enm.sp1$presence.points,
                                                       radius = 1000000, n = 10000, mask = climate_variables[[1]])
### M. arnholdi
enm.arnh <- enmtools.species()
enm.arnh$species.name <- "Marnholdi"
enm.arnh$presence.points <- arnh
enm.arnh$range <- background.raster.buffer(enm.arnh$presence.points, 50000, mask = climate_variables)
enm.arnh$background.points <- background.points.buffer(points = enm.arnh$presence.points,
                                                      radius = 1000000, n = 10000, mask = climate_variables[[1]])

### M. lehilahytsara
lehi <- na.omit(lehi)
enm.lehi <- enmtools.species()
enm.lehi$species.name <- "MlehilahytsaraS"
enm.lehi$presence.points <- lehi
enm.lehi$range <- background.raster.buffer(enm.lehi$presence.points, 50000, mask = climate_variables)
enm.lehi$background.points <- background.points.buffer(points = enm.lehi$presence.points,
                                                       radius = 1000000, n = 10000, mask = climate_variables[[1]])
### M. mittermeieri
enm.mitt <- enmtools.species()
enm.mitt$species.name <- "Mmittermeieri"
enm.mitt$presence.points <- mitt
enm.mitt$range <- background.raster.buffer(enm.mitt$presence.points, 50000, mask = climate_variables)
enm.mitt$background.points <- background.points.buffer(points = enm.mitt$presence.points,
                                                       radius = 1000000, n = 10000, mask = climate_variables[[1]])

### M. gerpi
enm.gerp <- enmtools.species()
enm.gerp$species.name <- "Mgerpi"
enm.gerp$presence.points <- gerp
enm.gerp$range <- background.raster.buffer(enm.gerp$presence.points, 50000, mask = climate_variables)
enm.gerp$background.points <- background.points.buffer(points = enm.gerp$presence.points,
                                                       radius = 1000000, n = 10000, mask = climate_variables[[1]])
### M. jollyae
enm.joll <- enmtools.species()
enm.joll$species.name <- "Mjollyae"
enm.joll$presence.points <- joll
enm.joll$range <- background.raster.buffer(enm.joll$presence.points, 50000, mask = climate_variables)
enm.joll$background.points <- background.points.buffer(points = enm.joll$presence.points,
                                                       radius = 1000000, n = 10000, mask = climate_variables[[1]])

### M. murinus N
enm.muriN <- enmtools.species()
enm.muriN$species.name <- "MmurinusN"
enm.muriN$presence.points <- muriN
enm.muriN$range <- background.raster.buffer(enm.muriN$presence.points, 50000, mask = climate_variables)
enm.muriN$background.points <- background.points.buffer(points = enm.muriN$presence.points,
                                                       radius = 1000000, n = 10000, mask = climate_variables[[1]])

### M. murinus C
enm.muriC <- enmtools.species()
enm.muriC$species.name <- "MmurinusC"
enm.muriC$presence.points <- muriC
enm.muriC$range <- background.raster.buffer(enm.muriC$presence.points, 50000, mask = climate_variables)
enm.muriC$background.points <- background.points.buffer(points = enm.muriC$presence.points,
                                                        radius = 1000000, n = 10000, mask = climate_variables[[1]])
# M. murinus south
enm.muriS <- enmtools.species()
enm.muriS$species.name <- "MmurinusS"
enm.muriS$presence.points <- muriS
enm.muriS$range <- background.raster.buffer(enm.muriS$presence.points, 50000, mask = climate_variables)
enm.muriS$background.points <- background.points.buffer(points = enm.muriS$presence.points,
                                                        radius = 1000000, n = 10000, mask = climate_variables[[1]])
### M. griseorufus
enm.gris <- enmtools.species()
enm.gris$species.name <- "Mgriseorufus"
enm.gris$presence.points <- gris
enm.gris$range <- background.raster.buffer(enm.gris$presence.points, 50000, mask = climate_variables)
enm.gris$background.points <- background.points.buffer(points = enm.gris$presence.points,
                                                       radius = 1000000, n = 10000, mask = climate_variables[[1]])


### M. margotmarshae
enm.marg <- enmtools.species()
enm.marg$species.name <- "Mmargotmarshae"
enm.marg$presence.points <- marg
enm.marg$range <- background.raster.buffer(enm.marg$presence.points, 50000, mask = climate_variables)
enm.marg$background.points <- background.points.buffer(points = enm.marg$presence.points,
                                                       radius = 1000000, n = 10000, mask = climate_variables[[1]])
### M. mamiratra
enm.mami <- enmtools.species()
enm.mami$species.name <- "Mmamiratra"
enm.mami$presence.points <- mami
enm.mami$range <- background.raster.buffer(enm.mami$presence.points, 50000, mask = climate_variables)
enm.mami$background.points <- background.points.buffer(points = enm.mami$presence.points,
                                                       radius = 1000000, n = 10000, mask = climate_variables[[1]])
### M. sambiranensis
enm.samb <- enmtools.species()
enm.samb$species.name <- "Msambiranensis"
enm.samb$presence.points <- samb
enm.samb$range <- background.raster.buffer(enm.samb$presence.points, 50000, mask = climate_variables)
enm.samb$background.points <- background.points.buffer(points = enm.samb$presence.points,
                                                       radius = 1000000, n = 10000, mask = climate_variables[[1]])

### M. tanosi
enm.tano <- enmtools.species()
enm.tano$species.name <- "Mtanosi"
enm.tano$presence.points <- tano
enm.tano$range <- background.raster.buffer(enm.tano$presence.points, 50000, mask = climate_variables)
enm.tano$background.points <- background.points.buffer(points = enm.tano$presence.points,
                                                       radius = 1000000, n = 10000, mask = climate_variables[[1]])


### M. tavaratra, northern and central populations
reduced_df1 <- dplyr::select(data, original_cluster, long, lat) # downsize dataframe
reduced_df1 <- na.omit(reduced_df1)

# reduce to species subsets
tavaN <- vect(filter(reduced_df1, original_cluster == "tavaN"), geom = c("long", "lat"), crs = "EPSG:4326")
tavaC <- vect(filter(reduced_df1, original_cluster == "tavaC"), geom = c("long", "lat"), crs = "EPSG:4326")


# tava north
enm.tavaN <- enmtools.species()
enm.tavaN$species.name <- "tavaN"
enm.tavaN$presence.points <- tavaN
enm.tavaN$range <- background.raster.buffer(enm.tavaN$presence.points, 50000, mask = climate_variables)
enm.tavaN$background.points <- background.points.buffer(points = enm.tavaN$presence.points,
                                                            radius = 1000000, n = 10000, mask = climate_variables[[1]])
# tava central
enm.tavaC <- enmtools.species()
enm.tavaC$species.name <- "tavaC"
enm.tavaC$presence.points <- tavaC
enm.tavaC$range <- background.raster.buffer(enm.tavaC$presence.points, 50000, mask = climate_variables)
enm.tavaC$background.points <- background.points.buffer(points = enm.tavaC$presence.points,
                                                        radius = 1000000, n = 10000, mask = climate_variables[[1]])





#### ENM estimation ####
# => the model fine tuning (stored in "args.xxx") object must be done beforehand.
#    See section "Model fine tuning" for this.
args.enm.rave = c("betamultiplier=4.0", "linear=TRUE", "quadratic=TRUE",
                  "product=FALSE", "hinge=FALSE", "threshold=FALSE", "autofeature=FALSE") # maxent parameters
rave.mx <- enmtools.maxent(enm.rave, climate_variables, test.prop = 0.0, args = args.enm.rave, clamp = FALSE) # maxent model

args.enm.danf = c("betamultiplier=2", "linear=TRUE", "quadratic=FALSE",
                  "product=TRUE", "hinge=FALSE", "threshold=FALSE", "autofeature=FALSE") # set maxent parameters
danf.mx <- enmtools.maxent(enm.danf, climate_variables, test.prop = 0.0, args = args.enm.danf, clamp = FALSE) # maxent modeling

args.enm.bong = c("betamultiplier=1", "linear=TRUE", "quadratic=FALSE",
                  "product=FALSE", "hinge=FALSE", "threshold=FALSE", "autofeature=FALSE") # maxent parameters
bong.mx <- enmtools.maxent(enm.bong, climate_variables, test.prop = 0.0, args = args.enm.bong, clamp = FALSE) # maxent model

args.enm.myox = c("betamultiplier=1", "linear=TRUE", "quadratic=FALSE",
                  "product=FALSE", "hinge=FALSE", "threshold=FALSE", "autofeature=FALSE") # maxent parameters
myox.mx <- enmtools.maxent(enm.myox, climate_variables, test.prop = 0.0, args = args.enm.myox, clamp = FALSE) # maxent model

args.enm.rufu = c("betamultiplier=1", "linear=TRUE", "quadratic=TRUE",
                  "product=FALSE", "hinge=FALSE", "threshold=FALSE", "autofeature=FALSE") # maxent parameters
rufu.mx <- enmtools.maxent(enm.rufu, climate_variables, test.prop = 0.0, args = args.enm.rufu, clamp = FALSE) # maxent model

args.enm.bert = c("betamultiplier=1.0", "linear=TRUE", "quadratic=TRUE",
                  "product=FALSE", "hinge=FALSE", "threshold=FALSE", "autofeature=FALSE") # maxent parameters
bert.mx <- enmtools.maxent(enm.bert, climate_variables, test.prop = 0.0, args = args.enm.bert, clamp = FALSE) # maxent model

args.enm.simm = c("betamultiplier=1", "linear=TRUE", "quadratic=FALSE",
                  "product=FALSE", "hinge=FALSE", "threshold=FALSE", "autofeature=FALSE") # maxent parameters
simm.mx <- enmtools.maxent(enm.simm, climate_variables, test.prop = 0.0, args = args.enm.simm, clamp = FALSE) # maxent model

args.enm.maca = c("betamultiplier=1", "linear=TRUE", "quadratic=TRUE",
                  "product=FALSE", "hinge=FALSE", "threshold=FALSE", "autofeature=FALSE") # maxent parameters
maca.mx <- enmtools.maxent(enm.maca, climate_variables, test.prop = 0.0, args = args.enm.maca, clamp = FALSE) # maxent model

args.enm.jona = c("betamultiplier=1", "linear=TRUE", "quadratic=TRUE",
                  "product=FALSE", "hinge=FALSE", "threshold=FALSE", "autofeature=FALSE") # maxent parameters
jona.mx <- enmtools.maxent(enm.jona, climate_variables, test.prop = 0.0, args = args.enm.jona, clamp = FALSE) # maxent model

args.enm.tava = c("betamultiplier=1", "linear=TRUE", "quadratic=FALSE",
                  "product=TRUE", "hinge=FALSE", "threshold=FALSE", "autofeature=FALSE") # maxent parameters
tava.mx <- enmtools.maxent(enm.tava, climate_variables, test.prop = 0.0, args = args.enm.tava, clamp = FALSE) # maxent model

args.enm.sp1 = c("betamultiplier=1", "linear=TRUE", "quadratic=TRUE",
                 "product=FALSE", "hinge=FALSE", "threshold=FALSE", "autofeature=FALSE") # maxent parameters
sp1.mx <- enmtools.maxent(enm.sp1, climate_variables, test.prop = 0.0, args = args.enm.sp1, clamp = FALSE) # maxent model

args.enm.arnh = c("betamultiplier=1", "linear=TRUE", "quadratic=FALSE",
                  "product=FALSE", "hinge=FALSE", "threshold=FALSE", "autofeature=FALSE") # maxent parameters
arnh.mx <- enmtools.maxent(enm.arnh, climate_variables, test.prop = 0.0, args = args.enm.arnh, clamp = FALSE) # maxent model

args.enm.lehi = c("betamultiplier=1", "linear=TRUE", "quadratic=TRUE",
                  "product=FALSE", "hinge=FALSE", "threshold=FALSE", "autofeature=FALSE") # maxent parameters
lehi.mx <- enmtools.maxent(enm.lehi, climate_variables, test.prop = 0.0, args = args.enm.lehi, clamp = FALSE) # maxent model

args.enm.mitt = c("betamultiplier=2.0", "linear=TRUE", "quadratic=TRUE",
                  "product=FALSE", "hinge=FALSE", "threshold=FALSE", "autofeature=FALSE") # maxent parameters
mitt.mx <- enmtools.maxent(enm.mitt, climate_variables, test.prop = 0.0, args = args.enm.mitt, clamp = FALSE) # maxent model

args.enm.gerp = c("betamultiplier=1", "linear=TRUE", "quadratic=FALSE",
                  "product=FALSE", "hinge=FALSE", "threshold=FALSE", "autofeature=FALSE") # maxent parameters
gerp.mx <- enmtools.maxent(enm.gerp, climate_variables, test.prop = 0.0, args = args.enm.gerp, clamp = FALSE) # maxent model

args.enm.joll = c("betamultiplier=1", "linear=TRUE", "quadratic=TRUE",
                  "product=FALSE", "hinge=FALSE", "threshold=FALSE", "autofeature=FALSE") # maxent parameters
joll.mx <- enmtools.maxent(enm.joll, climate_variables, test.prop = 0.0, args = args.enm.joll, clamp = FALSE) # maxent model

args.enm.muriN = c("betamultiplier=1", "linear=TRUE", "quadratic=FALSE",
                   "product=FALSE", "hinge=FALSE", "threshold=FALSE", "autofeature=FALSE") # maxent parameters
muriN.mx <- enmtools.maxent(enm.muriN, climate_variables, test.prop = 0.0, args = args.enm.muriN, clamp = FALSE) # maxent model

args.enm.muriC = c("betamultiplier=1", "linear=TRUE", "quadratic=TRUE",
                   "product=FALSE", "hinge=FALSE", "threshold=FALSE", "autofeature=FALSE") # maxent parameters
muriC.mx <- enmtools.maxent(enm.muriC, climate_variables, test.prop = 0.0, args = args.enm.muriC, clamp = FALSE) # maxent model

args.enm.mami = c("betamultiplier=2", "linear=TRUE", "quadratic=TRUE",
                  "product=FALSE", "hinge=FALSE", "threshold=FALSE", "autofeature=FALSE") # maxent parameters
mami.mx <- enmtools.maxent(enm.mami, climate_variables, test.prop = 0.0, args = args.enm.mami, clamp = FALSE) # maxent model

args.enm.marg = c("betamultiplier=2", "linear=TRUE", "quadratic=TRUE",
                  "product=FALSE", "hinge=FALSE", "threshold=FALSE", "autofeature=FALSE") # maxent parameters
marg.mx <- enmtools.maxent(enm.marg, climate_variables, test.prop = 0.0, args = args.enm.marg, clamp = FALSE) # maxent model

args.enm.samb = c("betamultiplier=1", "linear=TRUE", "quadratic=TRUE",
                  "product=FALSE", "hinge=FALSE", "threshold=FALSE", "autofeature=FALSE") # maxent parameters
samb.mx <- enmtools.maxent(enm.samb, climate_variables, test.prop = 0.0, args = args.enm.samb, clamp = FALSE) # maxent model

args.enm.tano = c("betamultiplier=1", "linear=TRUE", "quadratic=TRUE",
                  "product=FALSE", "hinge=FALSE", "threshold=FALSE", "autofeature=FALSE") # maxent parameters
tano.mx <- enmtools.maxent(enm.tano, climate_variables, test.prop = 0.0, args = args.enm.tano, clamp = FALSE) # maxent model

args.enm.gris = c("betamultiplier=1", "linear=TRUE", "quadratic=FALSE",
                  "product=TRUE", "hinge=FALSE", "threshold=FALSE", "autofeature=FALSE") # maxent parameters
gris.mx <- enmtools.maxent(enm.gris, climate_variables, test.prop = 0.0, args = args.enm.gris, clamp = FALSE) # maxent model

args.enm.muriS = c("betamultiplier=1", "linear=TRUE", "quadratic=FALSE",
                  "product=FALSE", "hinge=FALSE", "threshold=FALSE", "autofeature=FALSE") # maxent parameters
muriS.mx <- enmtools.maxent(enm.muriS, climate_variables, test.prop = 0.0, args = args.enm.muriS, clamp = FALSE) # maxent model

args.enm.tavaN = c("betamultiplier=1", "linear=TRUE", "quadratic=FALSE",
                   "product=FALSE", "hinge=FALSE", "threshold=FALSE", "autofeature=FALSE") # maxent parameters
tavaN.mx <- enmtools.maxent(enm.tavaN, climate_variables, test.prop = 0.0, args = args.enm.tavaN, clamp = FALSE) # maxent model

args.enm.tavaC = c("betamultiplier=1", "linear=TRUE", "quadratic=FALSE",
                   "product=TRUE", "hinge=FALSE", "threshold=FALSE", "autofeature=FALSE") # maxent parameters
tavaC.mx <- enmtools.maxent(enm.tavaC, climate_variables, test.prop = 0.0, args = args.enm.tavaC, clamp = FALSE) # maxent model




#### Niche overlap estimation ####
# empirical value with 95% confidence intervals.
# => to be calculated for each test case separately

# define test case; objects with ".mx" ending
test.tx1 <- rave.mx
test.tx2 <- danf.mx
niche_overlap <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(niche_overlap) <- c("Iteration", "env.D", "env.I")

for (i in 1:100) {
  overlap_result <- env.overlap(test.tx1, test.tx2, env = climate_variables, tolerance = 0.001,
                                max.reps = 99, cor.method = "spearman", chunk.size = 1e+02)
  overlap <- overlap_result[1:2]
  niche_overlap <- rbind(niche_overlap, as.numeric(c(i, overlap[1], overlap[2])))
}
colnames(niche_overlap) <- c("Iteration", "env.D", "env.I")


### calculate statistics for Schoener's D
mean(niche_overlap[,2])                  # mean value of 100 replicates
t.test(niche_overlap[,2])$conf.int[1:2]  # 95% confidence intervals







#### +++ Pre- and post-processing steps +++ ####

#### Model fine tuning and evaluation ####
# AIC-based parameter fine-tuning for maxent models; 
# do 1-by-1, just change data frame ("df_xxx") 
# extract AUC and CBI values of the best model for reporting model accuracy
validation <- ENMevaluate(df_rave[2:3], envs = climate_variables,
                          tune.args = list(fc = c("L", "P", "Q", "LP", "LQ", "H"),
                                           rm = 1:6),
                          n.bg = 10000, partitions = "randomkfold", # randomkfold with k = 5 (by default)
                          algorithm = 'maxnet', overlap = FALSE, 
                          bin.output = FALSE, doClamp = FALSE, 
                          parallel = TRUE, progbar = TRUE, numCores = 8)
validation <- as.data.frame(validation@results) # show results, choose model with deltaAIC = 0
validation # show results as df

write.xlsx(validation, file = "validation_rave.xlsx", sheetName="validation") # write results to table




#### Visual model evaluation ####
# => do 1-by-1
plot(rave.mx) # plotting for checking in geographic space
visualize.enm(rave.mx, climate_variables, 
              layers = c("PC1", "PC2"), plot.points = TRUE) # 2D plot for checking in environmental space
rave.mx$model # get variable contribution and plotting thresholds





#### Save ENMs as raster files ####
plot(rave.mx) # plot ENM to check validity
map.mx <- rave.mx$suitability # convert from enmtools-object to SpatRaster-object
writeRaster(map.mx, filename= "rave.mx", format="GTiff", overwrite = TRUE) # write as .tiff to working directory


