#! /bin/bash
#SBATCH -N 1
#SBATCH -n 2

source activate microc

bam=$1 #bam of Micro-C or MNase-seq data
coverage=$(echo $bam | sed 's/.bam/-coverage.txt/')

echo reporting total mapped reads
echo total mapped reads > $coverage
samtools view -c $bam >> $coverage

for mark in H3K9me3 H3K27me3 H3K27ac; do
	echo intersecting $mark
	bed=/ru-auth/local/home/iduba/linker-histone/Micro-C/annotations/chip/K562/K562-$mark\.st.bed
	intersected=$(basename "$bam" .bam)_$(basename "$bed" .bed).bam
	bedtools intersect -a $bam -b $bed -wa > $intersected
	echo reporting number reads for $mark
	echo $mark >> $coverage
	samtools view -c $intersected >> $coverage
done
 
