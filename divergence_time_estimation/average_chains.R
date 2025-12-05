#!/usr/bin/env Rscript
options(scipen = 999)  # Disable scientific notation

args <- commandArgs(trailingOnly = TRUE)
n <- length(args)

# Last argument is output file path
out <- args[n]

# All previous arguments are input chains
chain_files <- args[1:(n - 1)]

# Read the first file to initialize the sum
x_sum <- read.table(chain_files[1], header = TRUE, sep = "\t")

# Loop through remaining chains and add them
for (i in 2:length(chain_files)) {
  x <- read.table(chain_files[i], header = TRUE, sep = "\t")
  x_sum <- x_sum + x
}

# Average
x_avg <- x_sum / length(chain_files)

# Save to file
write.table(x_avg, file = out, row.names = FALSE, sep = "\t", quote = FALSE)
