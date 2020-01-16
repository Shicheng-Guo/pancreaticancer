
cd /gpfs/home/guosa/hpc/methylation/pancrease/medip/bam
mkdir temp

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/chromInfo.txt.gz -O chromInfo.txt.gz
wget https://raw.githubusercontent.com/Shicheng-Guo/Gscutility/master/bedGraph2wig.pl -O bedGraph2wig.pl
wget https://raw.githubusercontent.com/Shicheng-Guo/Gscutility/master/bdg2bw.sh -O bdg2bw.sh

for i in `ls *.bam`
do
echo $i
echo \#PBS -N $i  > $i.job
echo \#PBS -l nodes=1:ppn=6 >> $i.job
echo \#PBS -M Guo.shicheng\@marshfieldresearch.org >> $i.job
echo \#PBS -m abe  >> $i.job
echo \#PBS -o $(pwd)/temp/ >>$i.job
echo \#PBS -e $(pwd)/temp/ >>$i.job
echo cd $(pwd) >> $i.job
# echo bowtie2 -p 6 -x /gpfs/home/guosa/hpc/db/hg19/bowtie2/hg19 -1 $i\_1.fastq.gz -2 $i\_2.fastq.gz -S $i.sam >> $i.job
# echo samtools view -bS $i.sam \> $i.bam >> $i.job
# echo samtools sort $i.bam -o $i.sorted.bam -@6 >> $i.job
# echo samtools mpileup -uf ~/hpc/db/hg19/hg19.db $i.sorted.bam \| bcftools view -Ov - \> $i.bcf >> $i.job
# echo samtools depth $i.sorted.bam \> $i.wig >> $i.job
echo macs2 callpeak -t $i -f BAM -g hs -n $i -B -q 0.01 --SPMR >> $i.job
# echo  macs2 bdgcmp -t $i\_treat_pileup.bdg -c $i\_control_lambda.bdg -o $i.FE.bdg -m FE >>$i.job
# echo  macs2 bdgcmp -t $i\_treat_pileup.bdg -c $i\_control_lambda.bdg -o $i.logLR.bdg -m logLR -p 0.01 >>$i.job
# echo sh bdg2bw.sh $i.FE.bdg chromInfo.txt >>$i.job
# echo bedtools sort -i $i\_L003.sorted.macs_peaks.narrowPeak.bed \> $i.sort.bed >> $i.job
# echo bedGraphToBigWig $i.sort.bed chromInfo.txt $i.bw >>$i.job
# echo bigWigAverageOverBed $i.bw driven.hg19.bed $i.tab >>$i.job
# echo bigWigToWig ./ucsc/$i.FE.bw ./ucsc/$i.wig >> $i.job
# echo perl bedGraph2wig.pl --bedgraph $i.bam.FE.bdg --wig $i.wig --step 150 >> $i.job
# echo gzip -c $i.wig \> $i.wig.gz >> $i.job
# echo samtools sort PAAD.N.bam PAAD.N.sort.bam >> $i.job
qsub  $i.job
done

