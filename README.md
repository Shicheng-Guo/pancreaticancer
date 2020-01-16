## Genetic and Epigenetic of Pancreatic adenocarcinoma


```
cd ~/hpc/methylation/pancrease/medip/bam
```

Timeline: 

* 2020/01/13: renew bigwig and save in `/gpfs/home/guosa/hpc/methylation/pancrease/medip/bam`
* 2020/01/12: Update MBD-seq gene enhancer/promoter gene region methylation status
* 2020/01/11: TCGA gene expression analysis
* 2020/01/10: TCGA mutation analysis
* 2020/01/10: TCGA methylation analysis
* 2020/01/10: TCGA miRNA analysis


Method: 

* Macs2 was then used to generate read count normalized genome wide pileup tracks and lambda tracks for precipitated samples and input.(callpeak --nomodel --extsize 150 --SPMR) Pileup tracks were then Input corrected using the Macs2 subtract function (bdgcmp -m subtract). Log10 Fold enrichment was then calculated for each sample comparing input corrected pileups to input lambda. (bdgcmp -m logFE). Peaks were then called as regions with greater than 3 fold enrichment over lambda using the Macs2 bdgpeakcall function. (bdgpeakcall -c 0.477 -l 100 -g 100)


