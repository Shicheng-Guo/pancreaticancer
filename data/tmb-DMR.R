# cd ~/hpc/methylation/pancrease/medip
# cp ../pancrease.merge.sort.bed ./
# for i in 2019032901 2019032903 2019040901 2019051703 2019052301 2019053101 2019053102
# do
# cat $i\_*.venn >> pancrease.temp.bed
# done
# bedtools sort -i pancrease.temp.bed > pancrease.sort.bed
# bedtools merge -d 500 -i pancrease.sort.bed > pancrease.merge.bed
# rm pancrease.temp.bed
# rm pancrease.sort.bed
# awk '{print $1,$2,$3,$1":"$2"-"$3}' OFS="\t" pancrease.merge.bed > pancrease.merge.sort.bed
# cd ~/hpc/methylation/pancrease/medip/matrix
# for i in `ls *.bw`
# do
# bigWigAverageOverBed $i pancrease.merge.sort.bed $i.full.tab 
# echo $i.full.tab 
# done

file=list.files(pattern="*.tab")
data<-c()
for(i in 1:length(file)){
  temp<-read.table(file[i])
  data<-cbind(data,temp[,6])
  print(i)
}
colnames(data)<-unlist(lapply(strsplit(file,"_201907"), function(x) x[1]))
rownames(data)<-temp[,1]
file=c("2019032901","2019032903","2019040901","2019051703","2019052301","2019053101","2019053102")
DD<-list()
for (i in 1:length(file)){
  print(i)
  temp<-data[,grep(file[i],colnames(data))]+1
  ND<-temp[,grep("_N",colnames(temp))]
  CD<-temp[,-grep("_N",colnames(temp))]
  head(ND)
  head(CD)
  XX<-CD/ND
  DD[[i]]=XX
}
DD<-data.frame(DD,check.names=F)
tmb2<-read.table("https://raw.githubusercontent.com/Shicheng-Guo/PancreaticCancer/master/data/TMB_2.txt",head=T)
tmby<-tmb[match(colnames(DD),tmb2$TumoT_Sample_BaTcode),4]
PP<-c()
for(i in 1:nrow(DD)){
  fit=summary(lm(tmby~as.numeric(DD[i,])))
  if(nrow(fit$coefficients)>1){
    P<-fit$coefficients[2,]
  }else{
    P<-c(NA,NA,NA,NA)  
  }
  PP<-rbind(PP,P)
  print(i)
}
head(PP)
out=data.frame(as.character(cufflinks[,1]),PP,FDR=p.adjust(PP[,4],"fdr"))
head(out)
write.table(out,file="TMB-DMR.out.txt",sep="\t",quote = F,col.names = T,row.names = F)
library("qqman")
library("Haplin")
pdf("qqplot.pdf")
pvalues=na.omit(out[,5])
pQQ(pvalues, nlabs =length(pvalues), conf = 0.95)
dev.off()





