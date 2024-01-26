#! /bin/bash
#SBATCH -N 1
#SBATCH -n 1

source activate encode-atac-seq-pipeline
#tss=~/linker-histone/Micro-C/annotations/HG38_TSS_intersectby_K562ATACpeak_ConservativeIDR.bed
ctcf=/rugpfs/fs0/risc_lab/scratch/hcanaj/K562_vplot_byATACdata_061523/hg38_CTCF_byK562ChIPseq_bedformat.bed

for bw in *rmdup.bw; do
#	mat=`echo $bw | sed 's/.bw/-ATACTSS.mat.gz/'`
#	profile=`echo $bw | sed 's/.bw/-ATACTSSpileup.pdf/'`
	mat=`echo $bw | sed 's/.bw/-CTCF.mat.gz/'`
	profile=`echo $bw | sed 's/.bw/-CTCFpileup.pdf/'`
	computeMatrix reference-point -S $bw -R $ctcf -o $mat
	plotProfile -m $mat -o $profile --refPointLabel CTCF 
done
