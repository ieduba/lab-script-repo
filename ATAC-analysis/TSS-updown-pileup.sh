#! /bin/bash
#SBATCH -N 1
#SBATCH -n 1

source activate encode-atac-seq-pipeline

for bw in *bw; do
	matname=$bw-updownTSS.mat.gz
	profilename=$bw-updownTSS-pileup.pdf
	before=500
	after=1500
	computeMatrix reference-point -S $bw -R  WTvslowup-TSS-lfcsort.bed WTvslowdown-basemeanmatchnochange-TSS-lfcsort.bed \
	-o $matname -b $before -a $after
	plotProfile -m $matname -o $profilename --refPointLabel TSS 
done
