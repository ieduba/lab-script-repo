#! /bin/bash
#SBATCH -N 1
#SBATCH -n 1

source activate encode-atac-seq-pipeline

for mark in low scr; do
	for region in TSS activeTSS CTCF; do
		bws=`ls $mark\_*_in_human.bw`
		matname=$mark-$region.mat.gz
		profilename=$mark-$region\pileup.pdf
		
		if [ $region == TSS ]; then
			bed=~/linker-histone/Micro-C/annotations/hg38_TSS_export_bedrearrange_gene.bed
			before=500
			after=1500
		elif [ $region == activeTSS ]; then
			bed=~/linker-histone/Micro-C/annotations/HG38_TSS_intersectby_K562ATACpeak_ConservativeIDR.bed
			before=500
			after=1500
		else 
			bed=/rugpfs/fs0/risc_lab/scratch/hcanaj/K562_vplot_byATACdata_061523/hg38_CTCF_byK562ChIPseq_bedformat.bed
			before=1000
			after=1000
		fi

		computeMatrix reference-point -S $bws -R $bed -o $matname -b $before -a $after
		plotProfile -m $matname -o $profilename --refPointLabel $region --perGroup
	done
done
