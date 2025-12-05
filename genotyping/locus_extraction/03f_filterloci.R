#!/usr/bin/env Rscript

#### Script adapted and modified from Poelstra et al. 2021, Systematic Biology (https://doi.org/10.1093/sysbio/syaa053)

cat("\n#### 03f_filterloci.R:  Starting script.\n")

## Command-line arguments
options(warn = 1)
args <- commandArgs(trailingOnly = TRUE)

infile_locusstats <- args[1]
indir_fasta <- args[2]
outdir_fasta <- args[3]
maxmiss <- as.double(args[4])
mindist <- as.double(args[5])

cat("#### 03f_filterloci.R:  Input file with locus statistics:", infile_locusstats, "\n")
cat("#### 03f_filterloci.R:  Maximum percentage of missing data (maxmiss):", maxmiss, "\n")
cat("#### 03f_filterloci.R:  Minimum distance (bp) between loci (mindist):", mindist, "\n")
cat("#### 03f_filterloci.R:  Fasta input directory:", indir_fasta, "\n")
cat("#### 03f_filterloci.R:  Fasta output directory:", outdir_fasta, "\n")

## Other variables
if (!dir.exists(outdir_fasta)) dir.create(outdir_fasta, recursive = TRUE)

nrfiles <- length(list.files(indir_fasta))
cat("\n#### 03f_filterloci.R:  Number of files in indir_fasta:", nrfiles, "\n")

## Packages
library(pacman)
packages <- c("valr", "tidyverse")
p_load(char = packages, install = TRUE)


## Process input files
# Locus statistics
lstats <- read.delim(infile_locusstats, as.is = TRUE) %>% select(1:12)
colnames(lstats) <- c(
    "locus.full", "nInd", "bp", "nCells", "nN", "pN",
    "nvar", "pvar", "nPars", "pPars", "AT", "GC"
)
lstats <- lstats %>%
    mutate(locus = gsub(".fa", "", locus.full)) %>%
    separate(locus, into = c("scaffold", "pos"), sep = ":") %>%
    separate(pos, into = c("start", "end"), sep = "-") %>%
    mutate(start = as.integer(start), end = as.integer(end)) %>%
    select(
        locus.full, scaffold, start, end, nInd, bp, nN, pN, nvar, pvar, nPars,
        pPars, AT, GC
    ) %>%
    arrange(scaffold, start) %>%
    group_by(scaffold) %>%
    mutate(
        distToNext = lead(start) - (start + bp), # Include dist-to-next locus
        called_sites = round((bp * nInd) * (100 - pN))
    )

lstats_bed <- lstats %>%
    dplyr::select(scaffold, start, end, locus.full) %>%
    dplyr::rename(chrom = scaffold)

cat("\n#### 03f_filterloci.R:  Quantiles of locus length:\n")
print(quantile(lstats$bp))


## Filter for missing data
# Number of loci with certain amount of missing data
nrow.filter <- function(threshold) {
    lstats %>%
        filter(pN < threshold) %>%
        nrow()
}
cat("#### 03f_filterloci.R:  Number of loci:", nrow(lstats), "\n")
cat("#### 03f_filterloci.R:  Number of loci with <10% N:", nrow.filter(10), "\n")
cat("#### 03f_filterloci.R:  Number of loci with <5% N:", nrow.filter(5), "\n")
cat("#### 03f_filterloci.R:  Number of loci with <1% N:", nrow.filter(1), "\n")
cat("#### 03f_filterloci.R:  Number of loci with no Ns:", nrow.filter(0.001), "\n")

cat("\n#### 03f_filterloci.R:  Quantiles of missing data:\n")
quantile(lstats$pN)

miss_hi_rm <- lstats$locus.full[lstats$pN > maxmiss]


## Filter for close proximity
tooclose_locus1_index <- which(lstats$distToNext < mindist)
tooclose_locus2_index <- tooclose_locus1_index + 1
tooclose <- cbind(
    lstats[tooclose_locus1_index, c("locus.full", "called_sites")],
    lstats[tooclose_locus2_index, c("locus.full", "called_sites")]
)
colnames(tooclose) <- c("locus1", "called_sites1", "locus2", "called_sites2")

tooclose_rm <- unique(c(
    tooclose$locus2[which(tooclose$called_sites1 >= tooclose$called_sites2)],
    tooclose$locus1[which(tooclose$called_sites2 > tooclose$called_sites1)]
))

cat(
    "\n#### 03f_filterloci.R:  Number of loci to remove due to close proximity:",
    length(tooclose_rm), "\n"
)


## Copy good loci
badloci <- unique(c(tooclose_rm, miss_hi_rm))
cat("\n#### 03f_filterloci.R:  Total number of loci to remove:", length(badloci), "\n")

if (length(badloci) > 0) {
    loci_ok <- lstats %>%
        filter(!locus.full %in% badloci) %>%
        pull(locus.full)
}
if (length(badloci) == 0) loci_ok <- lstats %>% pull(locus.full)
cat("\n#### 03f_filterloci.R:  Number of loci selected:", length(loci_ok), "\n")
cat("\n#### 03f_filterloci.R:  First 10 loci:\n")
print(head(loci_ok))

loci_ok_files <- paste0(indir_fasta, "/", loci_ok)
nfiles_found <- sum(file.exists(loci_ok_files))
cat("\n#### 03f_filterloci.R:  Number of files found:", nfiles_found, "\n")

cat("\n#### 03f_filterloci.R:  Copying files to final directory ...\n")
loci_copied_files <- paste0(outdir_fasta, "/", loci_ok)
file.copy(from = loci_ok_files, to = loci_copied_files, overwrite = TRUE)


## Report:
cat("\n#### 03f_filterloci.R:  Counting number of files ...\n")
nrfiles <- length(list.files(outdir_fasta))
cat("\n#### 03f_filterloci.R:  Number of files in indir_fasta:", nrfiles, "\n")

cat("\n\n#### 03f_filterloci.R:  Done with script.\n")