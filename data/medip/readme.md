# medip-seq analysis in @CHG1

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
* DMR located in driver gene regions (PANCAN and PAAD)
* DMR located in promoter of driver gene (PANCAN and PAAD)
```
wget https://raw.githubusercontent.com/Shicheng-Guo/PancreaticCancer/master/data/medip/DriverMutationGene2018.txt

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


