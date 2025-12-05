#############################################
#Script for Acoustic data
#written by Marina Scheumann, 19.12.2022
#R-version: 4.2.2; Rstudio- version:2022.07.1 Build 554
#published in van Elst et al. (2024): Integrative taxonomy clarifies the evolution of a cryptic primate clade
#############################################

#########################
#LOAD packages
#########################
library(xlsx)           # v. 0.6.5 - Data import
library(dplyr)          # v. 1.1.4 - Data re-arrangement, summarizing, grouping
library(ggplot2)        # v. 3.5.1 - Plotting
library(ggpubr)         # v. 0.6.0 - plot
library(data.table)     # v. 1.15.4 - Summarize data
library(stats)          # v. 4.2.2 - PCA
library(factoextra)     # v. 1.07. - PCA eigenvalues
library(dynRB)          # v. 0.18 - Overlap analysis
library(nlme)           # v. 3.1-160 - LME Equipment
#library(devtools)      # necessary to install from github - Option 2: Cran packages only
#devtools::install_github("cmartin/ggConvexHull")
library(ggConvexHull)   # Polygon in plot


################################################################################
###########################
#Grey mouse lemur group
###########################
###########################
#IMPORT, PREPARE and CHECK DATA FILE
###########################
setwd("E:/00_Marina Festplatte/02_Marina Arbeit 2015/04_Projekte/19_Ute/05_Supplementary") # set new working directory") # set new working directory
df <- xlsx::read.xlsx("E:/00_Marina Festplatte/02_Marina Arbeit 2015/04_Projekte/19_Ute/05_Supplementary/AcousticData.xlsx",  
                            sheetName = "Murinus group_Advertisment")
View(df)
str(df)

    ####################################
    #### FOR SAFETY - Define and Check Factors
    ####################################
df$Location <- as.factor(df$Location) # define factor
class(df$Location)# check factor
levels(df$Location)  # check available location 
str(df)#temporal parameter in [ms]; frequency parameter in [kHz]

########################################
#Test for the effect of different recording equipment
#######################################
#LME for Microcebus (north) - comparison between analog and digital
#######################################
data_Amp <- subset(df, Location =="Ampijoroa")
View(data_Amp)

model1<-lme (Total_Duration ~ Recording_System, random=~1|Individual, data=data_Amp, method='ML')
summary (model1)

model1<-lme (Maxf0 ~ Recording_System, random=~1|Individual, data=data_Amp, method='ML')
summary (model1)

######################################
#Plot Equipment
######################################
#Order of the Location in the legend
data_E <- df %>% mutate(Location=recode(Location, 'Ampijoroa'='1', 'Kirindy'='2', 'Hannover'='3'))
data_E$Location <- as.character (data_E$Location)
class(data_E$Location)
View(data_E)

data_E$Total_Duration_ms <- data_E$Total_Duration*1000
data_E$Maxf0_kHz <- data_E$Maxf0/1000

totdur <- ggplot(data_E, aes(Location, Total_Duration_ms, fill=Recording_System))+
  geom_boxplot()+theme_bw()+
  scale_x_discrete(name = "Species", labels = c("M. murinus\n (north)", "M. murinus\n (central)", "M. ganzhorni"))+
  scale_y_continuous(name="Total duration of the Trill call [ms]")+
  theme(legend.spacing.y = unit(0.4, "cm"))+
  theme_bw()+
  theme(legend.position = "bottom",
        legend.title = element_text(size=12),
        axis.text = element_text(size=12),
        axis.title = element_text(size=14))+
  labs(fill="Recording system")
totdur

totmaxf0 <- ggplot(data_E, aes(Location, Maxf0_kHz, fill=Recording_System))+
  geom_boxplot()+theme_bw()+
  scale_x_discrete(name = "Species", labels = c("M. murinus\n(north)", "M. murinus\n(central)", "M. ganzhorni"))+
  scale_y_continuous(name="Maximum F0 of the Trill call [kHz]")+
  theme(legend.spacing.y = unit(0.4, "cm"))+
  theme_bw()+
  theme(legend.position = "bottom",
        legend.title = element_text(size=12),
        axis.text = element_text(size=12),
        axis.title = element_text(size=14))+
  labs(fill="Recording system")
totmaxf0

ggarrange (totdur, totmaxf0, labels="auto") 

ggsave("Equipment_murinus.tiff", plot = last_plot(), width = 12, height = 7, units = "in",
       dpi = 300, device = "tiff", scale = 1) # save plot as tiff file

##########################################
#PCA ANALYSIS
##########################################
    #Preparing data frame for analysis
    # reduce data to columns for PCA
    ####################################
pca.data <- dplyr::select(df,   # dataframe
                          Location,   # grouping variable (e.g., taxon labels)
                          Total_Duration, Duration, Meanf0,Sdf0, Minf0, Maxf0, VOI,Timemin, Timemax, Start_Meanf0, Bandwidth, # explanatory variables
                          Bandwidt_start.max, Slope)
pca.data <- na.omit(pca.data)  # drop NA values
table(pca.data$Ort)# see whats left in the dataframe
View(pca.data) # check dataframe
str(pca.data)
summary(pca.data)

    ##################################
    #perform PCA analysis
    ##################################
pca <- prcomp(pca.data[,2:length(pca.data)], # from 2 to length of pca.data 
              center = TRUE,
              scale. = TRUE)
eig.val <- get_eigenvalue(pca)
eig.val
summary(pca)# show PCA results

        ###################################
        #For safety - save acoustic parameters with PCA scores
        #COMBINE PCA FACTORS WITH GROUPING VARIABLE (Location) and acoustic parameters
        ####################################
data.P <- as.data.frame(cbind(pca.data, pca$x)) 
View(data.P)
summary (data.P)

write.csv(data.P
          ,file = "Grey_PCA_R_all calls.csv"
          ,quote = FALSE
          ,row.names = FALSE
)

################################
#overlap analysis with PCA scores
################################

      #########################
      #Prepare Table with Location and PCA Scores
      #########################

data.PO <- as.data.frame(cbind(df$Location, pca$x)) # combine it to define V1
View(data.PO)
summary (data.PO)

#PCA Daten
r <- dynRB_VPa(data.PO, steps = 51, pca.corr = FALSE, correlogram = TRUE)# max of 51 steps (range boxes) to reach conclusion
res_P <- r$result
res_P

##########
#Final Plot for summary sheets
############
#Order of the Location in the legend
data_N <- data.P %>% mutate(Location=recode(Location, 'Ampijoroa'='1', 'Kirindy'='2', 'Hannover'='3'))
data_N$Location <- as.character (data_N$Location)
class(data_N$Location)
View(data_N)

#Label and colour of Location
mylabels <- c(expression(paste(italic("M. murinus (north)"))),
              expression(paste(italic("M. murinus (central)"))),
              expression(paste(italic("M. ganzhorni"))))

mycolors <- c("#E69F00", "#56B4E9", "#009E73") # give colors, N-S direction of taxa

#Plotting; only change taxon labels ("Location") in aes for ggplot and aes for geom_convexhull
ggplot(data = data_N, aes(x=as.numeric(PC1), y=as.numeric(PC2), color=as.factor(Location))) + # Vx must be number of length
  geom_point(size=4, shape=20, show.legend=FALSE) +
  scale_color_manual(values=mycolors, name=" ", guide=guide_legend(ncol=3), labels=mylabels) +
  theme(panel.background = element_rect(fill = "white", colour = "black"), legend.position="bottom") +
  geom_convexhull(alpha = 0.05, aes(color = as.factor(Location))) +
  xlab("PC1 (32.6%)") + ylab("PC2 (24.4%)")

ggsave("Grey_PCA_End.tiff", plot = last_plot(), width = 6.35, height = 8.67, units = "in",
       dpi = 300, device = "tiff", scale = 1) # save plot as tiff file

################################################################################
###########################
#Ravelobenis-bongolavensis-danfossi
###########################
###########################
#IMPORT, PREPARE and CHECK DATA FILE
###########################
df2 <- xlsx::read.xlsx("E:/00_Marina Festplatte/02_Marina Arbeit 2015/04_Projekte/19_Ute/05_Supplementary/AcousticData.xlsx",  
                      sheetName = "Ravelo_bongo_danf_alert")
View(df2)
str(df2) #temporal parameter in [ms]; frequency parameter in [kHz]

####################################
#### FOR SAFETY - Define and Check Factors
####################################
df2$Species <- as.factor(df2$Species) # define factor
class(df2$Species)# check factor
levels(df2$Species)  # check available location 

##########################################
#PCA ANALYSIS
##########################################
#Preparing data frame for analysis
#reduce data to columns for PCA
####################################
pca.data2 <- dplyr::select(df2,   # dataframe
                          Species,   # grouping variable (e.g., taxon labels)
                          Duration, VOI, Minf0, Maxf0, Bandwidth, MeanF0, Sdf0, meanslope)
pca.data2 <- na.omit(pca.data2)  # drop NA values
table(pca.data2$Species)# see whats left in the dataframe
View(pca.data2) # check dataframe
str(pca.data2)
summary(pca.data2)

##################################
#perform PCA analysis
##################################
pca2 <- prcomp(pca.data2[,2:length(pca.data2)], # from 2 to length of pca.data 
              center = TRUE,
              scale. = TRUE)
eig.val <- get_eigenvalue(pca2)
eig.val
summary(pca2)# show PCA results

###################################
#For safety - save acoustic parameters with PCA scores
#COMBINE PCA FACTORS WITH GROUPING VARIABLE (Species) and acoustic parameters
####################################
data.P2 <- as.data.frame(cbind(pca.data2, pca2$x)) 
View(data.P2)
summary (data.P2)

write.csv(data.P2
          ,file = "RavBonDan_PCA_R_all calls.csv"
          ,quote = FALSE
          ,row.names = FALSE
)

################################
#overlap analysis with PCA scores
################################

#########################
#Prepare Table with Species and PCA Scores
#########################

data.PO2 <- as.data.frame(cbind(df2$Species, pca2$x)) # combine it to define V1
View(data.PO2)
summary (data.PO2)

#PCA data
r2 <- dynRB_VPa(data.PO2, steps = 51, pca.corr = FALSE, correlogram = TRUE)# max of 51 steps (range boxes) to reach conclusion
res_P2 <- r2$result
res_P2

##########
#Final Plot for summary sheets
############
#Order of the Species in the legend
data_N2 <- data.P2 %>% mutate(Species=recode(Species, 'danfossi'='1', 'bongolavensis'='2', 'ravelobensis'='3'))
data_N2$Species <- as.character (data_N2$Species)
class(data_N2$Species)
View(data_N2)

#Label and colour of Location
mylabels <- c(expression(paste(italic("M. danfossi"))),
              expression(paste(italic("M. bongolavensis"))),
              expression(paste(italic("M. ravelobensis"))))

mycolors <- c("#E69F00", "#56B4E9", "#009E73")# give colors, N-S direction of taxa

#Plotting; only change taxon labels ("Species") in aes for ggplot and aes for geom_convexhull
ggplot(data = data_N2, aes(x=as.numeric(PC1), y=as.numeric(PC2), color=as.factor(Species))) + # Vx must be number of length
  geom_point(size=4, shape=20, show.legend=FALSE) +
  scale_color_manual(values=mycolors, name=" ", guide=guide_legend(ncol=3), labels=mylabels) +
  theme(panel.background = element_rect(fill = "white", colour = "black"), legend.position="bottom") +
  geom_convexhull(alpha = 0.05, aes(color = as.factor(Species))) +
  xlab("PC1 (50.75%)") + ylab("PC2 (25.04%)")

ggsave("RavBonDan_PCA_End.tiff", plot = last_plot(), width = 6.35, height = 8.67, units = "in",
       dpi = 300, device = "tiff", scale = 1) # save plot as tiff file
