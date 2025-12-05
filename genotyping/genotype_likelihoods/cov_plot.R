#!/usr/bin/env Rscript

#### Script adapted and modified from J. Salmona

cat("\n#### cov_plot.R: Starting script.\n")

## Command-line arguments
options(warn = 1)
args <- commandArgs(trailingOnly = TRUE)
inds_file <- args[1]
bam_hits <- args[2]
set_id <- args[3]

## Report
cat("\n#### cov_plot.R: File with individuals:", inds_file, "\n")
cat("#### cov_plot.R: Directory with BAM hits:", bam_hits, "\n")
cat("#### cov_plot.R: Set ID:", set_id, "\n\n")

## Process command-line arguments
inds <- read.table(paste0(inds_file),header=F)
setwd(bam_hits)

## Overwrite output
system("mkdir -p ./statistics")
system(paste0("rm ./statistics/",set_id,".minmax.txt; touch ./statistics/",set_id,".minmax.txt"))

## Estimate and plot distribution and quantiles
cat("#### cov_plot.R: Estimating and plotting distribution ... \n")
pdf(paste0("./statistics/",set_id,".sbf1.f1.cov.pdf"),15,15)
par(mfrow=c(3,3))

for(indv in inds[,1]){
cat("#### cov_plot.R: Reading individual", indv, "\n")

tab <- read.table(paste0(indv,".sbf1.f1.bamhits"),header=F)
tab[,1] <- as.numeric(tab[,1])
tabx <- tab[tab[,1]>1,]
medcovx <- median(tabx[,1]) ;
tab2 <- tabx[tabx[,1]<medcovx*6,]
maxi <- round(quantile(tab2[,1], probs=0.996),digits=1)
mini <- round(quantile(tab2[,1], probs=0.005),digits=1)
medcov <- median(tab2[,1])

system(paste0("echo ",indv," ", mini," ", maxi," >> ./statistics/",set_id,".minmax.txt"))

hist1<-hist(tab2[,1],xlim=c(0,50),ylim=c(0,1e+4),breaks=seq(1,max(tab2[,1]),1),col="green",xlab="Coverage",main=indv)

  text(maxi+1,3.5e+3,paste("max cov",maxi),col="blue",srt=90)
  text(15,9e+3,paste(sum(hist1$counts[mini:maxi]),"loci with cov >",mini,"x and <",maxi,"x"),col="darkgreen")
  text(15,1e+4,paste("median cov =",medcov),col="red")
  abline(v=maxi,col="blue")
  text(0,5e+3,paste(length(tab[tab[,1]<2,1]),"loci with cov=1 (not shown)"),col="blue",srt=90)
}
dev.off()

## Report:
cat("\n#### cov_plot.R: Done with script.\n")