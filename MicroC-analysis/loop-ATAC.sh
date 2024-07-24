#! /bin/bash
#SBATCH -N 1
#SBATCH -n 4

source activate encode-atac-seq-pipeline

for ATAC in scr-KA23 dH1-KA23; do
	for loop in scronly lowonly; do
		bw=$ATAC\.bw
		matname=$ATAC-$loop.mat.gz
		profilename=$ATAC-$loop-tornado.pdf
		
		if [ $loop == scronly ]; then
			bed=~/linker-histone/Micro-C/in-situ/old-new-combo/sep-bio-reps/loops/diffloops-ATAC/scr-loop-anchors.bed
		else 
			bed=~/linker-histone/Micro-C/in-situ/old-new-combo/sep-bio-reps/loops/diffloops-ATAC/dH1-loop-anchors.bed
		fi
		computeMatrix reference-point -S $bw -R $bed -o $matname --referencePoint center -b 2500 -a 2500
		plotHeatmap -m $matname -o $profilename --refPointLabel loop
	done
done
