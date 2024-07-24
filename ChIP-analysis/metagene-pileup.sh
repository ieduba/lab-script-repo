#! /bin/bash
#SBATCH -N 1
#SBATCH -n 1

source activate encode-atac-seq-pipeline

for mark in low scr; do
	bws=`ls $mark\_*_in_human.bw`
	matname=$mark-metagene.mat.gz
	profilename=$mark-metagene-pileup.pdf

	computeMatrix scale-regions -S $bws -R ~/linker-histone/Micro-C/annotations/hg38_genebody_bedrearrange_gene.bed -o $matname -b 1000 -a 1000
	plotProfile -m $matname -o $profilename --perGroup 
done
