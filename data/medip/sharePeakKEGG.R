## loading packages
#BiocManager::install("ChIPseeker")
#BiocManager::install("seqplots")
#BiocManager::install("genomation")
#BiocManager::install("clusterProfiler")

setwd("~/hpc/methylation/pancrease/medip")
library("ChIPseeker")
library("seqplots")
library("genomation")
library("clusterProfiler")
library("ReactomePA")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

## KEGG enrichment analysis to shared DMRs (promoter and overall)
setwd("~/hpc/methylation/pancrease/medip/Intervene_results/sets")
REF1<-"1111_2019032903_T1.venn_2019032903_T2.venn_2019032903_T3.venn_2019032903_T4.venn.bed.ref"
REF2<-"1111_2019051703_T1.venn_2019051703_T2.venn_2019051703_T3.venn_2019051703_T4.venn.bed.ref"
REF3<-"1111_2019052301_T1.venn_2019052301_T2.venn_2019052301_T3.venn_2019052301_T4.venn.bed.ref"
REF4<-"1111_2019053101_T1.venn_2019053101_T2.venn_2019053101_T3.venn_2019053101_T4.venn.bed.ref"
REF5<-"1111_2019053102_T1.venn_2019053102_T2.venn_2019053102_T3.venn_2019053102_T4.venn.bed.ref"
REF6<-"111_2019032901_T1.venn_2019032901_T2.venn_2019032901_T3.venn.bed.ref"
REF7<-"111_2019040901_T1.venn_2019040901_T2.venn_2019040901_T3.venn.bed.ref"
for(i in c(REF1,REF2,REF3,REF4,REF5,REF6,REF7)){
  print(i)
  data<-read.table(i,head=F,sep="\t")
  set1<-subset(data,V10=="Promoter")
  write.table(set1,file=paste(i,"promoter.bed",sep=""),quote=F,sep="\t",col.names = F,row.names = F)
  peak <- readPeakFile(paste(i,"promoter.bed",sep=""),head=F)
  gene <- seq2gene(peak, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
  pathway <- enrichPathway(gene=gene,pvalueCutoff=0.05, readable=T)
  write.table(as.data.frame(pathway),file=paste(i,".dmr.peak.kegg.txt",sep=""),sep="\t",quote=F)  
  try(barplot(pathway, showCategory=20))
  ggsave(paste(i,".shared.promoter.peak.kegg1.pdf",sep=""))
  try(emapplot(pathway))
  ggsave(paste(i,".shared.promoter.peak.kegg2.pdf",sep=""))
  
  # full
  peak <- readPeakFile(i,head=F)
  gene <- seq2gene(peak, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
  pathway <- enrichPathway(gene=gene,pvalueCutoff=0.05, readable=T)
  try(barplot(pathway, showCategory=20))
  ggsave(paste(i,".shared.full.peak.kegg1.pdf",sep=""))
  try(emapplot(pathway))
  ggsave(paste(i,".shared.full.peak.kegg2.pdf",sep=""))
}

