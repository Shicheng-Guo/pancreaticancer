# CHG1
setwd("~/hpc/methylation/pancrease/medip")
file=c("2019032901","2019032903","2019040901","2019051703","2019052301","2019053101","2019053102")
for (i in file){
print(i)
input<-read.table(paste(i,".matrix",sep=""),head=T,row.names=1,check.names=F)
colnames(input)<-unlist(lapply(strsplit(colnames(input),"_20190"),function(x) x[1]))
foldc<-unlist(apply(input,1,function(x) mean(x[2:4]/x[1])))
foldc<-na.omit(foldc)
foldc<-foldc[-which(!is.finite(foldc))]
xsel<-head(order(foldc,decreasing=T),n=2000)
pdf(paste(i,".hvar.matrix.pdf",sep=""))
temp<-input[match(names(foldc)[xsel],rownames(input)),]
z<-apply(temp,1,function(x) (x-mean(x))/sd(x))
HeatMap(data.matrix(t(z)))
dev.off()
}
