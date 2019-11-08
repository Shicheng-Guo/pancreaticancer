setwd("//mcrfnas2/bigdata/Genetic/Projects/shg047/methylation/Pancreatic_cancer")
tmb<-read.csv("https://raw.githubusercontent.com/Shicheng-Guo/PancreaticCancer/master/data/TMB-1106.csv")
cufflinks<-read.table("cufflinksOut.txt",head=T,sep="\t",check.names = F)
dim(cufflinks)
newcolnames<-unlist(lapply(strsplit(colnames(cufflinks),"FRAL"),function(x) x[1]))
colnames(cufflinks)<-newcolnames
head(tmb)
head(cufflinks)
tail(sort(table(cufflinks$Symbol)))
idv<-unique(unlist(lapply(strsplit(newcolnames,"-"),function(x) x[1])))[2:8]
mean(tmb[,4])
DD<-list()
for(i in 1:length(idv)){
  temp<-cufflinks[,c(grep(idv[i],colnames(cufflinks)))]
  ND<-temp[,grep("-N",colnames(temp))]
  CD<-temp[,-grep("-N",colnames(temp))]
  XX<-CD-ND
  DD[[i]]=XX
}
DD<-data.frame(DD,check.names=F)
head(DD)

tmby<-tmb[match(colnames(DD),tmb$TumoT_Sample_BaTcode),4]
fit<-apply(DD,1,function(x) )

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
out=data.frame(cufflinks[,1],PP)
head(out)
write.table(out,file="TMB-DEG.out.txt",sep="\t",quote = F,col.names = NA,row.names = T)

