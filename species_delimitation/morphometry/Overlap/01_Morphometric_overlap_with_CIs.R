### Script by Dominik Schuessler

# In this script, morphological overlaps in multidimensional space are calculated, 
# based on four commonly shared variables across the whole genus.
# Taxa are abbreviated with a four letter code, based in the first four letters 
# of the species epithet. 

# packages
library(dynRB)          # calculating hypervolume overlap
library(readxl)         # data import
library(dplyr)          # data re-arrrangement, summarizing, grouping



# load data
setwd("C://Users//domin//OneDrive//Dokumente//011-Revised_taxonomy_of_Microcebus//Morpho")
data <- read_excel("01_Microcebus_morphological_data.xlsx", sheet = "data", col_names = TRUE)


# prepare data set
str(data)
data$population <- as.factor(data$population)


# rename to fit names in phylogeny/nodes
levels(data$population)
data$population <- recode(data$population, gris_Beza = 'gris', 
                         lehi_S = 'lehi',
                         lehi_N = 'lehi',
                         lehi_mitt = 'mitt')

# extract all taxa to single dataframes, to look if downsampling is needed; combine after; e.g., bora, maro not possible due to low sample size
arno <- dplyr::filter(data, population == "arno")
bert <- dplyr::filter(data, population == "bert")
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
muri_C <- dplyr::filter(data, population == "muri_C")
muri_N <- dplyr::filter(data, population == "muri_N")
myox <- dplyr::filter(data, population == "myox")
rave <- dplyr::filter(data, population == "rave")
rufu <- dplyr::filter(data, population == "rufu")
samb <- dplyr::filter(data, population == "samb")
simm <- dplyr::filter(data, population == "simm")
sp1 <- dplyr::filter(data, population == "sp1")
tano <- dplyr::filter(data, population == "tano")
tava <- dplyr::filter(data, population == "tava")

muri_N <- sample_n(muri_N, 150) # re-sampling to reduce sample size
rave <- sample_n(rave, 150) 
tava <- sample_n(tava, 100) 

# combine all taxa to one dataframe
data1 <- rbind(arno, bert, ganz, gerp, gris, bong, danf, joll, jona, mitt, lehi, maca,
               mami, marg, muri_C, muri_N, myox, rave, rufu, samb, simm, sp1, tano, tava) # combine without data deficient taxa (bora, maro)

 

#### Permutation tests of hypervolumes ####
# => define taxa to compare first, then run script and extract values at the end
data2 <- select(data1, population, ear.length, head.length,
                body.mass, tail.length)     # select the variables of interest
data3 <- na.omit(data2)
data4 <- droplevels(data3)
str(data4)
data4$population <- as.factor(data4$population)
table(data4$population)



# define the test case (species pair)
tx1 <- dplyr::filter(data4, population == "danf") # taxon 1 to test
tx2 <- dplyr::filter(data4, population == "rave") # taxon 2 to test

results_list <- vector("list", 99)

# Loop 99 times
for (i in 1:99) {
  # For tx1 dataframe
  num_rows_to_extract_tx1 <- round(nrow(tx1) * 9/10) # 10-fold approach, 1/10 is always left out for calculation
  selected_rows_tx1 <- sample(nrow(tx1), num_rows_to_extract_tx1)
  extracted_dataframe_tx1 <- tx1[selected_rows_tx1, ]
  
  # For tx2 dataframe
  num_rows_to_extract_tx2 <- round(nrow(tx2) * 9/10) # 10-fold approach, 1/10 is always left out for calculation
  selected_rows_tx2 <- sample(nrow(tx2), num_rows_to_extract_tx2)
  extracted_dataframe_tx2 <- tx2[selected_rows_tx2, ]
  
  # Combine the extracted dataframes
  test_group <- rbind(extracted_dataframe_tx1, extracted_dataframe_tx2)
  test_group <- droplevels(test_group)
  
  # Perform your analysis on 'test_group'
  overlap <- dynRB_VPa(test_group, 
                       steps = 51, pca.corr = FALSE, 
                       correlogram = FALSE)$result[2:3, 1:5]
  
  # Save the result in the results_list
  results_list[[i]] <- overlap[which.max(overlap$port_prod), ]
}

result_dataframe <- do.call(rbind, results_list)

hist(result_dataframe$port_prod, xlab = "Morphological overlap estimate", main = "") # visualize estimate distribution
CI <- quantile(result_dataframe$port_prod, c(0.05, 0.95)) # 95% confidence intervals


# empirical value
test_group <- droplevels(rbind(tx1, tx2))
empirical <- max(dynRB_VPa(test_group, 
                           steps = 51, pca.corr = FALSE, 
                           correlogram = FALSE)$result[2:3, 1:5]$port_prod)

## results
empirical # empirical value (maximized morphological overlap)
CI        # 95% confidence interval from jackknife resampling








#### reference case: M. tavaratra ####
tava <- droplevels(tava)

levels(as.factor(tava$Place))       # available places 
tava$Place <- recode(tava$Place,    # recode sampling locations to split in northern and central/southern populations
                     Ambilondambo = 'tavaC', 
                     Ambohitsitondroina = 'tavaC',
                     Ampanetibe_Madirobe = 'tavaN',
                     Ampasimaty = 'tavaC',
                     Analabe = 'tavaN',
                     Analahefina = 'tavaC',
                     Analamerana = 'tavaC',
                     Andavakoera_1 = 'tavaC',
                     Andavakoera_2 = 'tavaC',
                     Andrafiambany = 'tavaC',
                     Andranomadiro = 'tavaN',
                     Ankaramy = 'tavaC',
                     Ankarana = 'tavaC',
                     Ankavana = 'tavaN',
                     Antsahabe = 'tavaC',
                     Antsaharaingy = 'tavaC',
                     Antsakay = 'tavaC',
                     Bekaraoka = 'tavaC',
                     Benanofy = 'tavaC',
                     Binara = 'tavaC',
                     Bobankora = 'tavaC',
                     Bobankora_est = 'tavaC',
                     Menagisy = 'tavaC',
                     Orangea = 'tavaN',
                     Analafiana = 'tavaC',
                     Solaniampilana = 'tavaC') 
tava_split <- dplyr::select(tava, Place, ear.length, head.length,
                      body.mass, tail.length)
table(tava_split$Place)


tx1 <- dplyr::filter(tava_split, Place == "tavaN") # taxon 1 to test
tx2 <- dplyr::filter(tava_split, Place == "tavaC") # taxon 2 to test

results_list <- vector("list", 99)

# Loop 99 times
for (i in 1:99) {
  # For tx1 dataframe
  num_rows_to_extract_tx1 <- round(nrow(tx1) * 9/10) # 10-fold approach, 1/10 is always left out for calculation
  selected_rows_tx1 <- sample(nrow(tx1), num_rows_to_extract_tx1)
  extracted_dataframe_tx1 <- tx1[selected_rows_tx1, ]
  
  # For tx2 dataframe
  num_rows_to_extract_tx2 <- round(nrow(tx2) * 9/10) # 10-fold approach, 1/10 is always left out for calculation
  selected_rows_tx2 <- sample(nrow(tx2), num_rows_to_extract_tx2)
  extracted_dataframe_tx2 <- tx2[selected_rows_tx2, ]
  
  # Combine the extracted dataframes
  test_group <- rbind(extracted_dataframe_tx1, extracted_dataframe_tx2)
  test_group <- droplevels(test_group)
  
  # Perform your analysis on 'test_group'
  overlap <- dynRB_VPa(test_group, 
                       steps = 51, pca.corr = FALSE, 
                       correlogram = FALSE)$result[2:3, 1:5]
  
  # Save the result in the results_list
  results_list[[i]] <- overlap[which.max(overlap$port_prod), ]
}

result_dataframe <- do.call(rbind, results_list)

hist(result_dataframe$port_prod, xlab = "Morphological overlap estimate", main = "") # visualize estimate distribution
CI <- quantile(result_dataframe$port_prod, c(0.05, 0.95)) # 95% confidence intervals


# empirical value
test_group <- droplevels(rbind(tx1, tx2))
empirical <- max(dynRB_VPa(test_group, 
                           steps = 51, pca.corr = FALSE, 
                           correlogram = FALSE)$result[2:3, 1:5]$port_prod)

## results
empirical # empirical value (maximized morphological overlap)
CI        # 95% confidence interval from jackknife resampling





