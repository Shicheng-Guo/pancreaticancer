cd ~/hpc/methylation/pancrease/medip
wget https://raw.githubusercontent.com/Shicheng-Guo/PancreaticCancer/master/data/medip/DriverMutationGene2018.txt
###################################################################################################
## Venn for promoter and enhancer regions to PAAD and PANCAN driver mutation genes
grep PAAD DriverMutationGene2018.txt | awk '{print $1}' > paad.txt
grep PANCAN DriverMutationGene2018.txt | awk '{print $1}' > pancan.txt

grep -w -f paad.txt ~/hpc/db/hg19/refGeneV2.hg19.bed | awk '{print $1,$2,$3,$6,$7}' OFS="\t" | sort -u > paad.TumorDrivenMutationList.hg19.bed
grep -w -f pancan.txt ~/hpc/db/hg19/refGeneV2.hg19.bed | awk '{print $1,$2,$3,$6,$7}' OFS="\t" | sort -u > pancan.TumorDrivenMutationList.hg19.bed

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

