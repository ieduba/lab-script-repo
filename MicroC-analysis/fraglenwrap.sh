#! /bin/bash
#SBATCH -N 1
#SBATCH -n 2

source activate microc

bam=$1
sam=`basename $bam | sed 's/.bam/.sam/'`
lens=`basename $bam | sed 's/.bam/-len.txt/'`
hist=`basename $bam | sed 's/.bam/-lenhist.pdf/'`

echo making sam

samtools view $bam > $sam

python fraglen.py -s $sam -o $lens -p $hist
