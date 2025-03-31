#! /bin/bash
#SBATCH -N 1
#SBATCH -n 2

source activate data

bam=$1
bigwig=`echo $bam | sed 's/.bam/.bw/'`

bamCoverage -b $bam -o $bigwig --filterRNAstrand forward -p 2
