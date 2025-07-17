#!/usr/bin/env Rscript
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

# Parse arguments
wd <- args[1]
mcmc_file <- args[2]
out_file <- args[3]
mutrate <- as.numeric(args[4])
mutrate_var <- ifelse(length(args) >= 11, as.numeric(args[5]), NA)
gentime <- as.numeric(args[6])
gentime_sd <- ifelse(length(args) >= 11, as.numeric(args[7]), NA)
burnin <- as.integer(args[8])
last_sample <- as.integer(args[9])
distribution <- tolower(args[10]) %in% c("true", "t", "1")

# Set working directory
setwd(wd)

# Function to process MCMC log: remove burn-in and truncate if needed
cut_log <- function(logfile, burnin, last_sample) {
  Log <- read.table(logfile, header = TRUE)
  cat('\n## Processing:', logfile, "\nInitial rows:", nrow(Log), "\n")

  if (burnin > 0) {
    Log <- Log[Log$Gen >= burnin, ]
    cat("After burn-in removal:", nrow(Log), "rows\n")
  }

  if (!is.na(last_sample)) {
    Log <- Log[Log$Gen <= last_sample, ]
    cat("After truncating to last sample:", nrow(Log), "rows\n")
  }

  if (nrow(Log) < 1000) {
    warning("Less than 1000 rows left after filtering — skipping.")
    return(NULL)
  }

  return(Log)
}

# Read and filter MCMC log
mcmc <- cut_log(mcmc_file, burnin, last_sample)
if (is.null(mcmc)) quit(save = "no", status = 1)

# Convert to data.frame
mcmc <- as.data.frame(mcmc)

# Identify theta and tau columns
theta_cols <- grep("theta", colnames(mcmc))
tau_cols <- grep("tau", colnames(mcmc))

# Apply either distribution or fixed scaling
if (distribution) {
  ## Generate distributions
  n <- nrow(mcmc)

  # Generation time: lognormal distribution
  gentime_dist <- rlnorm(n, meanlog = log(gentime), sdlog = log(gentime_sd))

  # Mutation rate: gamma distribution, scaled down to per-site
  shape <- mutrate^2 / mutrate_var
  rate <- mutrate / mutrate_var
  mutrate_dist <- rgamma(n, shape = shape, rate = rate) * 1e-8

  # Scale theta and tau
  mcmc[, theta_cols] <- mcmc[, theta_cols] / (4 * mutrate_dist)
  mcmc[, tau_cols] <- mcmc[, tau_cols] * gentime_dist / mutrate_dist

} else {
  # Use fixed point estimates (mutation rate in per-site units)
  mutrate <- mutrate * 1e-8

  mcmc[, theta_cols] <- mcmc[, theta_cols] / (4 * mutrate)
  mcmc[, tau_cols] <- mcmc[, tau_cols] * gentime / mutrate / 1000  # output in ky
}

# Save output
write.table(mcmc, file = out_file, row.names = FALSE, sep = "\t", quote = FALSE)