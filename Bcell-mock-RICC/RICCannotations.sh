#! /bin/bash
# SBATCH -N 1
# SBATCH -n 10 

echo This is job $SLURM_JOB_ID
echo Test permission

source activate "ATACseq"

bam=$1

cd ~/linker-histone/Bcell

annos=`ls *trim.bed`

echo 'bam:' $bam
for anno in $annos; do
	echo 'annotation:' $anno
	bamshort=`echo $bam | cut -d'_' -f1`
	annoshort=`echo $anno | cut -d'.' -f1`
	out=~/linker-histone/Bcell/$annoshort-$bamshort'.bam'
	bedtools intersect -abam $bam -b $anno  > $out
done
