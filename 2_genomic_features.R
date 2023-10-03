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

#from this point on in your R session only chromosome 1 would be
#consulted when you call the various retrieval methods… If you need 
#to reset back to the original seqlevels (i.e. to the seqlevels stored 
#in the db), then set the seqlevels to seqlevels0(txdb)

seqlevels(txdb) <- seqlevels0(txdb)

#use seqlevels to set only chromosome 15 to be active 

seqlevels(txdb) <- "chr15"
seqlevels(txdb)

seqlevels(txdb)

#Retrieving data using the select method

# The TxDb objects inherit from AnnotationDb objects (just as
#the ChipDb and OrgDb objects do). One of the implications of this 
#relationship is that these object ought to be used in similar ways
#to each other. Therefore we have written supporting columns, keytypes,
#keys and select methods for TxDb objects.
#
# These methods can be a useful way of extracting data from a TxDb 
#object. And they are used in the same way that they would be used
#to extract information about a ChipDb or an OrgDb object. Here is 
#a simple example of how to find the UCSC transcript names that match 
#with a set of gene IDs.

keys <- c("100033416", "100033417", "100033420")
columns(txdb)
# [1] "CDSCHROM"   "CDSEND"     "CDSID"      "CDSNAME"    "CDSSTART"  
# [6] "CDSSTRAND"  "EXONCHROM"  "EXONEND"    "EXONID"     "EXONNAME"  
# [11] "EXONRANK"   "EXONSTART"  "EXONSTRAND" "GENEID"     "TXCHROM"   
# [16] "TXEND"      "TXID"       "TXNAME"     "TXSTART"    "TXSTRAND"  
# [21] "TXTYPE"   

keytypes(txdb)
# [1] "CDSID"    "CDSNAME"  "EXONID"   "EXONNAME" "GENEID"   "TXID"    
# [7] "TXNAME" 

select(txdb, keys = keys, columns = "TXNAME", keytype="GENEID")
#     GENEID     TXNAME
# 1 100033416 uc001yxl.4
# 2 100033417 uc001yxo.3
# 3 100033420 uc001yxr.3

columns (txdb)
# [1] "CDSCHROM"   "CDSEND"     "CDSID"      "CDSNAME"    "CDSSTART"  
# [6] "CDSSTRAND"  "EXONCHROM"  "EXONEND"    "EXONID"     "EXONNAME"  
# [11] "EXONRANK"   "EXONSTART"  "EXONSTRAND" "GENEID"     "TXCHROM"   
# [16] "TXEND"      "TXID"       "TXNAME"     "TXSTART"    "TXSTRAND"  
# [21] "TXTYPE

cols <- c("TXNAME", "TXSTRAND", "TXCHROM")
select(txdb, keys=keys, columns=cols, keytype="GENEID")

# Methods for returning GRanges objects
# 
# Retrieving data with select is useful, but sometimes it is more convenient to extract the result as GRanges objects. This is often the case when you are doing counting or specialized overlap operations downstream. For these use cases there is another family of methods available.
# 
# Perhaps the most common operations for a TxDb object is to retrieve the genomic coordinates or ranges for exons, transcripts or coding sequences. The functions transcripts, exons, and cds return the coordinate information as a GRanges object.
# 
# As an example, all transcripts present in a TxDb object can be obtained as follows:

GR <- transcripts (txdb)
GR[1:3]
# GRanges object with 3 ranges and 2 metadata columns:
#   seqnames            ranges strand |     tx_id
# <Rle>         <IRanges>  <Rle> | <integer>
#   [1]    chr15 20362688-20364420      + |     53552
# [2]    chr15 20487997-20496811      + |     53553
# [3]    chr15 20723929-20727150      + |     53554
# tx_name
# <character>
#   [1]  uc001yte.1
# [2]  uc001ytf.1
# [3]  uc001ytj.3
# -------
#   seqinfo: 1 sequence from hg19 genome

# The transcripts function returns a GRanges class 
#object. You can learn a lot more about the manipulation
#of these objects by reading the GenomicRanges
#introductory vignette. The show method for a GRanges 
#object will display the ranges, seqnames (a chromosome
#or a contig), and strand on the left side and then 
#present related metadata on the right side. At the 
#bottom, the seqlengths display all the possible
#seqnames along with the length of each sequence.

# 
# The strand function is used to obtain the strand information from the transcripts. The sum of the Lengths of the Rle object that strand returns is equal to the length of the GRanges object.

tx_strand <- strand(GR)
tx_strand

sum(runLength(tx_strand))
#3337

length(GR)

# In addition, the transcripts function can also be used to retrieve a subset of the transcripts available such as those on the
#\(+\)-strand of chromosome 1.

