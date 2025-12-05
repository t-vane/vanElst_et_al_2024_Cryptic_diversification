#!/usr/bin/env Rscript

# Load required packages
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(data.table, tidyverse)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
logfile <- args[1]
outfile <- args[2]
burnin <- as.double(args[3])
last_sample <- as.double(args[4])

# Define population hierarchy (i.e., parent nodes/tree topology)
poplist <- list(Marnholdi = 'arn_sp1',  Mberthae = 'ber_ruf', Mbongolavensis = 'bon_rav', Mboraha = 'sim_bor', Mdanfossi = 'dan_bonrav', Mganzhorni = 'gan_man',
				Mgerpi = 'ger_marjol', Mgriseorufus = 'gri_mur', Mmarohita = 'jol_mar', Mjollyae = 'jol_mar', Mjonahi = 'jon_mac', MlehilahytsaraN = 'leh_mit', MlehilahytsaraS = 'leh', Mmacarthurii = 'jon_mac',
				Mmamiratra = 'mar_mam', Mmanitatra = 'gan_man', Mmargotmarshae = 'mar_mam', Mmittermeieri = 'leh_mit', MmurinusN = 'mur', MmurinusC = 'mur_ganman', Mmyoxinus = 'myo_berruf', Mravelobensis = 'bon_rav', 
				Mrufus = 'ber_ruf', Msambiranensis = 'sam_marmam', Msimmonsi = 'sim_bor', Msp1 = 'arn_sp1', Mtanosi = 'tan_leh', Mtavaratra = 'tav_arnsam', Mirzazaza = 'root', 
				arn_sp1 = 'arnsp1_sammarmam', ber_ruf = 'myo_berruf', bon_rav = 'dan_bonrav', sim_bor = 'simbor_tavetc', dan_bonrav = 'danbonrav_mur', gan_man = 'mur_ganman', ger_marjol = 'germarjol_simetc', 
				gri_mur = 'danbonrav_mur', jol_mar = 'ger_marjol', jon_mac = 'jonmac_east', leh_mit = 'leh', leh = 'leh_myoberruf', mar_mam = 'sam_marmam', mur = 'danbonrav_mur', mur_ganman = 'mur', 
				myo_berruf = 'leh_myoberruf', sam_marmam = 'arnsp1_sammarmam', tan_leh = 'tavarnsam_tanlehmyoberruf', tav_arnsam = 'tavarnsam_tanlehmyoberruf', 
				arnsp1_sammarmam = 'tavarnsam_tanlehmyoberruf', simbor_tavetc = 'germarjol_simetc', danbonrav_mur = 'micro', germarjol_simetc = 'jonmac_east', jonmac_east = 'micro', 
				leh_myoberruf = 'tan_leh', tavarnsam_tanlehmyoberruf = 'simbor_tavetc',
				germarjol_simetc = 'jonmac_east', micro = 'root')

# Function to cut the logfile
cut_log <- function(logfile, burnin, last_sample) {
  Log <- read.table(logfile, header = TRUE)
  Log <- Log[Log$Sample >= burnin, ]
  if (!is.null(last_sample)) Log <- Log[Log$Sample <= last_sample, ]
  return(Log)
}

# Function to calculate gdi estimates
add_gdi <- function(Log, poplist) {
  get_parent <- function(p, lvl) {
    for (i in seq_len(lvl)) p <- poplist[[p]]
    return(p)
  }

  get_gdi <- function(p, lvl) {
    parent <- get_parent(p, lvl)
    theta <- Log[[paste0("theta_", p)]]
    tau <- Log[[paste0("tau_", parent)]]
    1 - exp((-2 * tau) / theta)
  }

  for (lvl in 1:4) {
    for (p in head(names(poplist), 28)) {
      parent <- tryCatch(get_parent(p, lvl), error = function(e) NULL)
      if (!is.null(parent) &&
          paste0("theta_", p) %in% names(Log) &&
          paste0("tau_", parent) %in% names(Log)) {
        colname <- paste0("gdi", ifelse(lvl == 1, "", lvl), "_", p)
        Log[[colname]] <- get_gdi(p, lvl)
      }
    }
  }

  return(Log)
}

# Run the process
Log <- cut_log(logfile, burnin, last_sample)
if (!is.null(Log)) {
  Log <- add_gdi(Log, poplist)
  
  # Select only "Gen" + gdi columns
  gdi_cols <- grep("^gdi", names(Log), value = TRUE)
  out <- Log[, c("Gen", gdi_cols), drop = FALSE]
  
  # Write result
  outdir <- dirname(outfile)
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  write.table(out, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
}
