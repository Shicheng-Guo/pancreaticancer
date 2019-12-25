# medip-seq analysis in @CHG1

* MBD-Seq data for 33 samples of pancreatic cancer were downloaded from [SRP222713](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP222713&o=acc_s%3Aa&s=SRR10150657,SRR10150659,SRR10150661,SRR10150662,SRR10150663,SRR10150664,SRR10150665,SRR10150666,SRR10150667,SRR10150668,SRR10150669,SRR10150670,SRR10150671,SRR10150672,SRR10150673,SRR10150674,SRR10150675,SRR10150676,SRR10150677,SRR10150678,SRR10150679,SRR10150680,SRR10150681,SRR10150683,SRR10150684,SRR10150685,SRR10150686,SRR10150687,SRR10150688,SRR10150689,SRR10150660,SRR10150658,SRR10150682)

```
cd ~/hpc/methylation/pancrease/medip
```
MACS2
```
for i in `ls *sorted`
do
macs2 callpeak -t $i -f BAM -g hs -n $i.macs -B -q 0.05  &
done
```

narrowPeaks2bigwig
```
for i in `ls *macs_peaks.narrowPeak`
do 
awk '{print $1,$2,$3,$5}' OFS="\t" $i | bedtools sort -i - > ./matrix/$i.sort.bdg
bedGraphToBigWig ./matrix/$i.sort.bdg ~/hpc/db/hg19/hg19.chrom.sizes ./matrix/$i.bw 
echo $i 
done
```

bedtools to substract normal peaks from cancer genome.
```
bedtools intersect -a 2019032901_T1_20190709N_CCAGTTCA_S73_L003.sorted.macs_peaks.narrowPeak -b 2019032901_N_20190709N_CCTCCTGA_S76_L003.sorted.macs_peaks.narrowPeak -v > 2019032901_T1.venn
bedtools intersect -a 2019032901_T2_20190709N_CCGAAGTA_S74_L003.sorted.macs_peaks.narrowPeak -b 2019032901_N_20190709N_CCTCCTGA_S76_L003.sorted.macs_peaks.narrowPeak -v > 2019032901_T2.venn
bedtools intersect -a 2019032901_T3_20190709N_CCGTGAGA_S75_L003.sorted.macs_peaks.narrowPeak -b 2019032901_N_20190709N_CCTCCTGA_S76_L003.sorted.macs_peaks.narrowPeak -v > 2019032901_T3.venn
bedtools intersect -a 2019032903_T1_20190709N_CGAACTTA_S77_L003.sorted.macs_peaks.narrowPeak -b 2019032903_N_20190709N_CTGAGCCA_S81_L003.sorted.macs_peaks.narrowPeak -v > 2019032903_T1.venn
bedtools intersect -a 2019032903_T2_20190709N_CGACTGGA_S78_L003.sorted.macs_peaks.narrowPeak -b 2019032903_N_20190709N_CTGAGCCA_S81_L003.sorted.macs_peaks.narrowPeak -v > 2019032903_T2.venn
bedtools intersect -a 2019032903_T3_20190709N_CGCATACA_S79_L003.sorted.macs_peaks.narrowPeak -b 2019032903_N_20190709N_CTGAGCCA_S81_L003.sorted.macs_peaks.narrowPeak -v > 2019032903_T3.venn
bedtools intersect -a 2019032903_T4_20190709N_CTCAATGA_S80_L003.sorted.macs_peaks.narrowPeak -b 2019032903_N_20190709N_CTGAGCCA_S81_L003.sorted.macs_peaks.narrowPeak -v > 2019032903_T4.venn
bedtools intersect -a 2019040901_T1_20190709N_CTGGCATA_S82_L003.sorted.macs_peaks.narrowPeak -b 2019040901_N_20190709N_GAGCTGAA_S85_L003.sorted.macs_peaks.narrowPeak -v > 2019040901_T1.venn
bedtools intersect -a 2019040901_T2_20190709N_GAATCTGA_S83_L003.sorted.macs_peaks.narrowPeak -b 2019040901_N_20190709N_GAGCTGAA_S85_L003.sorted.macs_peaks.narrowPeak -v > 2019040901_T2.venn
bedtools intersect -a 2019040901_T3_20190709N_CAAGACTA_S84_L003.sorted.macs_peaks.narrowPeak -b 2019040901_N_20190709N_GAGCTGAA_S85_L003.sorted.macs_peaks.narrowPeak -v > 2019040901_T3.venn
bedtools intersect -a 2019051703_T1_20190709N_GATAGACA_S86_L003.sorted.macs_peaks.narrowPeak -b 2019051703_N_20190709N_GCTCGGTA_S90_L003.sorted.macs_peaks.narrowPeak -v > 2019051703_T1.venn
bedtools intersect -a 2019051703_T2_20190709N_GCCACATA_S87_L003.sorted.macs_peaks.narrowPeak -b 2019051703_N_20190709N_GCTCGGTA_S90_L003.sorted.macs_peaks.narrowPeak -v > 2019051703_T2.venn
bedtools intersect -a 2019051703_T3_20190709N_GCGAGTAA_S88_L003.sorted.macs_peaks.narrowPeak -b 2019051703_N_20190709N_GCTCGGTA_S90_L003.sorted.macs_peaks.narrowPeak -v > 2019051703_T3.venn
bedtools intersect -a 2019051703_T4_20190709N_GCTAACGA_S89_L003.sorted.macs_peaks.narrowPeak -b 2019051703_N_20190709N_GCTCGGTA_S90_L003.sorted.macs_peaks.narrowPeak -v > 2019051703_T4.venn
bedtools intersect -a 2019052301_T1_20190709N_GGAGAACA_S91_L003.sorted.macs_peaks.narrowPeak -b 2019052301_N_20190709N_GTCTGTCA_S95_L003.sorted.macs_peaks.narrowPeak -v > 2019052301_T1.venn
bedtools intersect -a 2019052301_T2_20190709N_GGTGCGAA_S92_L003.sorted.macs_peaks.narrowPeak -b 2019052301_N_20190709N_GTCTGTCA_S95_L003.sorted.macs_peaks.narrowPeak -v > 2019052301_T2.venn
bedtools intersect -a 2019052301_T3_20190709N_GTACGCAA_S93_L003.sorted.macs_peaks.narrowPeak -b 2019052301_N_20190709N_GTCTGTCA_S95_L003.sorted.macs_peaks.narrowPeak -v > 2019052301_T3.venn
bedtools intersect -a 2019052301_T4_20190709N_GTCGTAGA_S94_L003.sorted.macs_peaks.narrowPeak -b 2019052301_N_20190709N_GTCTGTCA_S95_L003.sorted.macs_peaks.narrowPeak -v > 2019052301_T4.venn
bedtools intersect -a 2019053101_T1_20190709N_GTGTTCTA_S96_L003.sorted.macs_peaks.narrowPeak -b 2019053101_N_20190709N_TCTTCACA_S100_L003.sorted.macs_peaks.narrowPeak -v > 2019053101_T1.venn
bedtools intersect -a 2019053101_T2_20190709N_TAGGATGA_S97_L003.sorted.macs_peaks.narrowPeak -b 2019053101_N_20190709N_TCTTCACA_S100_L003.sorted.macs_peaks.narrowPeak -v > 2019053101_T2.venn
bedtools intersect -a 2019053101_T3_20190709N_TATCAGCA_S98_L003.sorted.macs_peaks.narrowPeak -b 2019053101_N_20190709N_TCTTCACA_S100_L003.sorted.macs_peaks.narrowPeak -v > 2019053101_T3.venn
bedtools intersect -a 2019053101_T4_20190709N_TCCGTCTA_S99_L003.sorted.macs_peaks.narrowPeak -b 2019053101_N_20190709N_TCTTCACA_S100_L003.sorted.macs_peaks.narrowPeak -v > 2019053101_T4.venn
bedtools intersect -a 2019053102_T1_20190709N_TGAAGAGA_S101_L003.sorted.macs_peaks.narrowPeak -b 2019053102_N_20190709N_TTCACGCA_S105_L003.sorted.macs_peaks.narrowPeak -v > 2019053102_T1.venn
bedtools intersect -a 2019053102_T2_20190709N_TGGAACAA_S102_L003.sorted.macs_peaks.narrowPeak -b 2019053102_N_20190709N_TTCACGCA_S105_L003.sorted.macs_peaks.narrowPeak -v > 2019053102_T2.venn
bedtools intersect -a 2019053102_T3_20190709N_TGGCTTCA_S103_L003.sorted.macs_peaks.narrowPeak -b 2019053102_N_20190709N_TTCACGCA_S105_L003.sorted.macs_peaks.narrowPeak -v > 2019053102_T3.venn
bedtools intersect -a 2019053102_T4_20190709N_TGGTGGTA_S104_L003.sorted.macs_peaks.narrowPeak -b 2019053102_N_20190709N_TTCACGCA_S105_L003.sorted.macs_peaks.narrowPeak -v > 2019053102_T4.venn
```
venn shown for sample specific MACS peaks
```
for i in `ls *venn`
do
awk '{print $1,$2,$3}' OFS="\t" $i > $i.bed
done
conda install -c bioconda intervene

cd ~/hpc/methylation/pancrease/medip
for i in 2019032901 2019032903 2019040901 2019051703 2019052301 2019053101 2019053102
do
mkdir $i
intervene venn -i /gpfs/home/guosa/hpc/methylation/pancrease/medip/venn/$i*.bed --project Intervene_results
intervene upset -i /gpfs/home/guosa/hpc/methylation/pancrease/medip/venn/$i*.bed --project Intervene_results
intervene pairwise  -i /gpfs/home/guosa/hpc/methylation/pancrease/medip/venn/$i*.bed --project Intervene_results
done
intervene upset -i /gpfs/home/guosa/hpc/methylation/pancrease/medip/venn/*.bed --project Intervene_results
intervene pairwise  -i /gpfs/home/guosa/hpc/methylation/pancrease/medip/venn/*.bed --project Intervene_results
```
peakheatmap
```
cd ~/hpc/methylation/pancrease/medip
for i in 2019032901 2019032903 2019040901 2019051703 2019052301 2019053101 2019053102
do
cat $i\_*.venn > $i.temp.bed
bedtools sort -i $i.temp.bed > $i.sort.bed
bedtools merge -d 500 -i $i.sort.bed > $i.merge.bed
rm $i.temp.bed
rm $i.sort.bed
awk '{print $1,$2,$3,$1":"$2"-"$3}' OFS="\t" $i.merge.bed > $i.merge.sort.bed
for j in `ls $i\_*.bw`
do
bigWigAverageOverBed $j $i.merge.sort.bed $j.tab 
done 
done
```
tab2matrix
```
cd ~/hpc/methylation/pancrease/medip
for i in 2019032901 2019032903 2019040901 2019051703 2019052301 2019053101 2019053102
do
perl tab2matrix.pl $i > $i.matrix
done
```

Here, we want to check the DMRs in Driver Mutation Genes ()
* DMR located in driver gene regions ([200 PANCAN and 11 PAAD](DriverMutationGene2018.GeneList.txt))
* DMR located in promoter of driver gene ([200 PANCAN and 11 PAAD](DriverMutationGene2018.GeneList.txt))
```
cd ~/hpc/methylation/pancrease/medip
mkdir ./paad/
for i in `ls -d *.venn`
do
bedtools intersect -wa -a $i -b paad.TumorDrivenMutationList.hg19.bed > ./paad/$i.pancan.venn
awk '{print $1,$2,$3}' OFS="\t" ./paad/$i.pancan.venn > ./paad/$i.pancan.bed
done
cd ./paad/
# conda install -c bioconda intervene
for i in 2019032901 2019032903 2019040901 2019051703 2019052301 2019053101 2019053102
do
intervene venn -i $i*.bed --project $i.pancan.fullgene
done

cd ~/hpc/methylation/pancrease/medip
mkdir ./pancan/
for i in `ls -d *.venn`
do
bedtools intersect -wa -a $i -b pancan.TumorDrivenMutationList.hg19.bed > ./pancan/$i.pancan.venn
awk '{print $1,$2,$3}' OFS="\t" ./pancan/$i.pancan.venn > ./pancan/$i.pancan.bed
done
cd ./pancan/
# conda install -c bioconda intervene
for i in 2019032901 2019032903 2019040901 2019051703 2019052301 2019053101 2019053102
do
intervene venn -i $i*.bed --project $i.pancan.fullgene
done

###################################################################################################
## Venn for promoter and enhancer regions to PAAD and PANCAN driver mutation genes

grep -w -f paad.txt ~/hpc/db/hg19/refGeneV2.hg19.bed | awk '{print $1,$2,$3,$6,$7}' OFS="\t" | grep "Promoter\|Enhancer" | sort -u > paad.promoter.TumorDrivenMutationList.hg19.bed
grep -w -f pancan.txt ~/hpc/db/hg19/refGeneV2.hg19.bed | awk '{print $1,$2,$3,$6,$7}' OFS="\t" | grep "Promoter\|Enhancer" | sort -u > pancan.promoter.TumorDrivenMutationList.hg19.bed
cd ~/hpc/methylation/pancrease/medip
rm -rf ./paad/
mkdir ./paad/
for i in `ls -d *.venn`
do
bedtools intersect -wa -a $i -b paad.promoter.TumorDrivenMutationList.hg19.bed > ./paad/$i.paad.venn
awk '{print $1,$2,$3}' OFS="\t" ./paad/$i.paad.venn > ./paad/$i.paad.bed
done
cd ./paad/
# conda install -c bioconda intervene
for i in 2019032901 2019032903 2019040901 2019051703 2019052301 2019053101 2019053102
do
intervene venn -i $i*.bed --project $i.paad.promoterEnhancaer
done

cd ~/hpc/methylation/pancrease/medip
rm -rf ./pancan/
mkdir ./pancan/
for i in `ls -d *.venn`
do
bedtools intersect -wa -a $i -b pancan.promoter.TumorDrivenMutationList.hg19.bed > ./pancan/$i.pancan.venn
awk '{print $1,$2,$3}' OFS="\t" ./pancan/$i.pancan.venn > ./pancan/$i.pancan.bed
done
cd ./pancan/
# conda install -c bioconda intervene
for i in 2019032901 2019032903 2019040901 2019051703 2019052301 2019053101 2019053102
do
intervene venn -i $i*.bed --project $i.pancan.promoterEnhancaer
done
```

Gene Set heatmap 
```
cd ~/hpc/methylation/pancrease/medip

grep PAAD DriverMutationGene2018.txt | awk '{print $1}' > paad.txt
grep PANCAN DriverMutationGene2018.txt | awk '{print $1}' > pancan.txt
grep -w -f paad.txt ~/hpc/db/hg19/refGeneV2.hg19.bed | awk '{print $1,$2,$3,$6,$7}' OFS="\t" | sort -u > paad.TumorDrivenMutationList.hg19.bed
grep -w -f pancan.txt ~/hpc/db/hg19/refGeneV2.hg19.bed | awk '{print $1,$2,$3,$6,$7}' OFS="\t" | sort -u > pancan.TumorDrivenMutationList.hg19.bed

for i in 2019032901 2019032903 2019040901 2019051703 2019052301 2019053101 2019053102
do
cat $i\_*.venn > $i.temp.bed
bedtools sort -i $i.temp.bed > $i.sort.bed
bedtools merge -d 500 -i $i.sort.bed > $i.merge.bed
rm $i.temp.bed
rm $i.sort.bed
awk '{print $1,$2,$3,$1":"$2"-"$3}' OFS="\t" $i.merge.bed > $i.merge.sort.bed
for j in `ls $i\_*.bw`
do
bigWigAverageOverBed $j $i.merge.sort.bed $j.tab 
done 
done

###############################################################################################
## heatmap for all the samples of each individual

cd ~/hpc/methylation/pancrease/medip/heatmap/pancan
wget https://raw.githubusercontent.com/Shicheng-Guo/PancreaticCancer/master/data/medip/HeatMap.R -O HeatMap.R 
source("./HeatMap.R")
file=c("2019032901","2019032903","2019040901","2019051703","2019052301","2019053101","2019053102")
for (i in file){
print(i)
input<-read.table(paste(i,".matrix",sep=""),head=T,row.names=1,check.names=F)
colnames(input)<-unlist(lapply(strsplit(colnames(input),"_"),function(x) x[2]))
foldc<-unlist(apply(input,1,function(x) mean(x[2:4]/x[1])))
foldc<-na.omit(foldc)
if(sum(!is.finite(foldc))>0){
foldc<-foldc[-which(!is.finite(foldc))]
}
xsel<-head(order(foldc,decreasing=T),n=2000)
pdf(paste(i,".panca.fullgene.matrix.pdf",sep=""))
temp<-input[match(names(foldc)[xsel],rownames(input)),]
z<-apply(temp,1,function(x) (x-mean(x))/sd(x))
HeatMap(data.matrix(t(z)))
dev.off()
}

###############################################################################################
## heatmap for all the samples of each individual
matrix<-list.files(pattern="*.matrix$")
data<-read.table(matrix[1],head=T,row.names=1,check.names=F)
colnames(data)<-unlist(lapply(strsplit(colnames(data),"_201"),function(x) x[1]))
for(i in 2:length(matrix)){
print(i)
temp<-read.table(matrix[i],head=T,row.names=1,check.names=F)
colnames(temp)<-unlist(lapply(strsplit(colnames(temp),"_201"),function(x) x[1]))
data<-data.frame(data,temp,check.names=F)
}
head(data)

pdf("heatmap.panca.pdf")
data<-data.matrix(data.frame(data))
heatmap(cor(data),Rowv=F,Colv=F)
dev.off()


```


