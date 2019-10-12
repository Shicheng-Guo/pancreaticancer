## loading packages
BiocManager::install("ChIPseeker")
BiocManager::install("seqplots")
BiocManager::install("genomation")
BiocManager::install("clusterProfiler")

library("ChIPseeker")
library("seqplots")
library("genomation")
library("clusterProfiler")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library("clusterProfiler")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
