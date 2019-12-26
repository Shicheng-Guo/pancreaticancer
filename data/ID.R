sam1<-read.csv("https://raw.githubusercontent.com/Shicheng-Guo/PancreaticCancer/master/data/medip/saminfo.csv",head=F)
sam2<-read.csv("https://raw.githubusercontent.com/Shicheng-Guo/PancreaticCancer/master/data/medip/fastq/SraRunTable-SRP222713.csv")
sam<-data.frame(sam1,sam2[match(sam1$V4,sam2$Library.Name),])
write.csv(sam,file="newsam-mbdseq.csv")
