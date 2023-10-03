#Genomic Features analysis for CSDE1 RIP 

#https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html

#Queenie Tsang
#October 3 2023

#clear environment variables
rm(list=ls())

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicFeatures")

library("GenomicFeatures")

#tutorial walkthrough
#https://bioconductor.org/packages/release/bioc/vignettes/GenomicFeatures/inst/doc/GenomicFeatures.html

#loading Transcript Data
samplefile <- system.file("extdata", "hg19_knownGene_sample.sqlite",
                          package="GenomicFeatures")
txdb <- loadDb(samplefile)
txdb

#or more commonly just load a TxDb annotation package 
install.packages("TxDb.Hsapiens.UCSC.hg19.knownGene")

BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene 
txdb

#loading the package like this creates a TxDb object and 
# by default the object will have the same name as the package 
#itself 

#Prefiltering data based on Chromosomes

#it's possible to filter the data that's returned from a TxDb object
#based on its chromosome. This can be a useful way to limit the things that are returned if you are only interested in studying a handful of chromosomes.

#To determine which chromosomes are currently active, use the seqlevels
#method. For example:
head(seqlevels(txdb))
#[1] "chr1" "chr2" "chr3" "chr4" "chr5" "chr6"

# Will tell you all the chromosomes that are active for the TxDb.Hsapiens.UCSC.hg19.knownGene TxDb object (by default it will be all of them).
# 
# If you then wanted to only set Chromosome 1 to be active you could do it like this:
  
seqlevels(txdb) <- "chr1"


