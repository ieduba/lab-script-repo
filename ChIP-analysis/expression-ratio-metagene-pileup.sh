#! /bin/bash
#SBATCH -N 1
#SBATCH -n 4

source activate encode-atac-seq-pipeline

for ratio in low_10_in low_20_in scr_10_in scr_20_in; do
	bws=`ls $ratio\_human.bw`
	echo $bws
	matname=$ratio-expression-metagene.mat.gz
	profilename=$ratio-expression-metagene-pileup.pdf

	q1=~/linker-histone/RNA-seq/human/DESeq2_results/K562_scr1_q1_genes-genebody.bed
	q4=~/linker-histone/RNA-seq/human/DESeq2_results/K562_scr1_q4_genes-genebody.bed

	echo making $matname
	computeMatrix scale-regions -S $bws -R $q1 $q4 -o $matname -b 1000 -a 1000
	echo plotting $profilename
	plotProfile -m $matname -o $profilename 
done
