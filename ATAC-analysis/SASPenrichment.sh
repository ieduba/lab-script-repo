#! /bin/bash
# SBATCH -N 1
# SBATCH -n 10 

echo This is job $SLURM_JOB_ID
echo Test permission

source activate "ATACseq"

file=$1

cd LSA2/hg19-pipelineRun

echo 'sample:' $file
fileshort=`echo $file | cut -d'-' -f1`
bam=$fileshort/*.rmdup.bam

out=~/scripts/LSA2/SASPgenes/$fileshort-SASP-hg19.bam
bedtools intersect -a $bam -b ~/scripts/LSA2/SASPgenes/SASPgenes-hg19.final.bed -wa > $out

