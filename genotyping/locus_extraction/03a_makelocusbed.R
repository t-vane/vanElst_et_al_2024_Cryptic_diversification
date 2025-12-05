#!/usr/bin/env Rscript

#### Script adapted and modified from Poelstra et al. 2021, Systematic Biology (https://doi.org/10.1093/sysbio/syaa053)

cat("\n#### 03a_makelocusbed.R: Starting script.\n")

## Command-line arguments
options(warn = 1)
args <- commandArgs(trailingOnly = TRUE)
set_id <- args[1]
infile_inds <- args[2]
indir_bed <- args[3]
outfile_bed <- args[4]

min_element_ovl <- as.double(args[5])
min_element_ovl_trim <- as.double(args[6])
min_locus_size <- as.double(args[7])
max_dist_within_ind <- as.double(args[8])
max_dist_between_ind <- as.double(args[9])
min_element_size <- as.double(args[10])

last_row <- strtoi(args[11])

## Packages
library(pacman)
library(valr)
library(IRanges)
packages <- c("data.table", "tidyverse")
p_load(char = packages, install = TRUE)

## Process command-line arguments
ids <- readLines(infile_inds)
bedfiles <- paste0(indir_bed, "/", ids, ".callable.bed")

min_element_ovl <- length(bedfiles) * min_element_ovl
min_element_ovl_trim <- length(bedfiles) * min_element_ovl_trim

## Report
cat("\n#### 03a_makelocusbed.R: Set ID:", set_id, "\n")
cat("#### 03a_makelocusbed.R: Directory with BED files:", indir_bed, "\n")
cat("#### 03a_makelocusbed.R: BED files:", bedfiles, "\n")
cat("#### 03a_makelocusbed.R: Number of BED files:", length(bedfiles), "\n")
cat("#### 03a_makelocusbed.R: BED output file:", outfile_bed, "\n\n")

cat("#### 03a_makelocusbed.R: Maximum distance within individuals:", max_dist_within_ind, "\n")
cat("#### 03a_makelocusbed.R: Maximum distance between individuals:", max_dist_between_ind, "\n")
cat("#### 03a_makelocusbed.R: Minimum element overlap for locus creation:", min_element_ovl, "\n")
cat("#### 03a_makelocusbed.R: Minimum element overlap for locus trimming:", min_element_ovl_trim, "\n")
cat("#### 03a_makelocusbed.R: Minimum element size:", min_element_size, "\n")
cat("#### 03a_makelocusbed.R: Minimum locus size:", min_locus_size, "\n")

cat("#### 03a_makelocusbed.R: Number of loci to process (all if 0):", last_row, "\n")

## Set functions
get_indloci <- function(bedfile, max_dist, last_row = 0) {
    cat("#### 03a_makelocusbed.R: get_indloci function: bedfile:", bedfile, "\n")

    bed <- fread(bedfile,
        header = FALSE,
        colClasses = c("character", "integer", "integer", "character"),
        col.names = c("chrom", "start", "end", "status")
    )
    bed$status <- NULL

    if (last_row == 0) {
        cat("#### 03a_makelocusbed.R: get_indloci(): Using all rows ...\n")
    } else {
        cat(
            "#### 03a_makelocusbed.R: get_indloci(): Selecting rows until row:",
            last_row, "\n"
        )
        bed <- bed[1:last_row, ]
    }

    bed <- as.tbl_interval(bed)
    bed <- bed_merge(bed, max_dist = max_dist)
    bed <- arrange(bed, start)
    return(bed)
}

collect_indloci <- function(bedfiles, last_row,
                            max_dist_within_ind, min_element_size) {
    bedlist <- lapply(bedfiles, get_indloci,
        max_dist = max_dist_within_ind, last_row = last_row
    )
    bed_by_ind <- do.call(rbind, bedlist)

    # This step is necessary for bed_merge to work
    bed_by_ind <- data.frame(
        chrom = bed_by_ind$chrom,
        start = bed_by_ind$start,
        end = bed_by_ind$end
    ) %>%
        as.tbl_interval() %>%
        arrange(start) %>%
        filter(end - start >= min_element_size)

    return(bed_by_ind)
}

trim_locus <- function(row.nr, locus_df, element_df, min_element_ovl_trim) {
    locus <- locus_df[row.nr, ]
    locus_id <- paste0(locus$chrom, "_", locus$start)

    locus_elements <- bed_intersect(bed_by_ind, locus)
    bed_cov_ir <- IRanges(
        start = locus_elements$start.x,
        end = locus_elements$end.x
    )
    cov <- IRanges::coverage(bed_cov_ir)

    ok <- which(cov@values > min_element_ovl_trim)

    if (length(ok >= 1)) {
        first_ok <- ok[1]
        if (first_ok > 1) {
            first_base <- sum(cov@lengths[1:first_ok])
        } else {
            first_base <- locus$start
        }

        last_ok <- ok[length(ok)]
        if (last_ok < length(cov@values)) {
            last_base <- sum(cov@lengths[1:last_ok])
        } else {
            last_base <- locus$end
        }

        locus_length <- locus$end - locus$start
        trimmed_start <- first_base - locus$start
        trimmed_end <- locus$end - last_base
        locus_length_final <- last_base - first_base

        cat(
            row.nr, locus_id, "/ length:", locus_length,
            "/ trimmed start:", trimmed_start, "/ trimmed end:", trimmed_end,
            "/ remaining bases:", locus_length_final, "\n"
        )

        locus_trimmed <- data.frame(
            chrom = locus$chrom,
            start = first_base, end = last_base
        )
        locus_trimmed_length <- locus_trimmed$end - locus_trimmed$start

        if (locus_trimmed_length < min_locus_size) {
            cat(row.nr, "Locus too small...\n")
        } else {
            return(locus_trimmed)
        }
    } else {
        cat(row.nr, "#### 03a_makelocusbed.R: Coverage too low: skipping locus ...\n")
    }
}

## Get per-individual bed file and merge loci
bed_by_ind <- collect_indloci(bedfiles,
    last_row = last_row,
    max_dist_within_ind = max_dist_within_ind,
    min_element_size = min_element_size
)

## Merge per-individual loci
bed_merged <- bed_merge(bed_by_ind, max_dist = max_dist_between_ind)
cat("\n## Nr of initial loci:", nrow(bed_merged), "\n")

## Calculate "coverage" (number of elements) overlapping with each locus
## Retain only those with "min.elements" number of overlapping elements, and "min.frac" fraction of overlap (latter is not very important)
bed_cov <- bed_coverage(bed_merged, bed_by_ind) %>%
    filter(.ints >= min_element_ovl)
cat(
    "#### 03a_makelocusbed.R: Number of loci after filtering by coverage:",
    nrow(bed_cov), "\n\n"
)

## Trim loci
cat("#### 03a_makelocusbed.R: Trimming loci ...\n")
bed_trim_list <- lapply(1:nrow(bed_cov), trim_locus,
    locus_df = bed_cov,
    element_df = bed_by_ind,
    min_element_ovl_trim = min_element_ovl_trim
)
bed_trim <- do.call(rbind, bed_trim_list)
cat("#### 03a_makelocusbed.R: Trimming done.\n\n")
cat("#### 03a_makelocusbed.R: Number of loci retained:", nrow(bed_trim), "\n\n")

## Save bedfile
cat("#### 03a_makelocusbed.R: Writing bedfile ...", outfile_bed, "\n")
write.table(bed_trim, outfile_bed,
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)

## Report:
cat("\n#### 03a_makelocusbed.R: Done with script.\n")