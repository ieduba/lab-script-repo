#! /bin/bash
# SBATCH -N 1
# SBATCH -n 2 

echo This is job $SLURM_JOB_ID
echo Test permission

source activate "ATACseq"

directory=$1
cd ~/linker-histone/MESR3/corrected
cd $directory
pwd

sample=`basename $directory | sed 's:/::'` 

merged=$sample-merge.bam
FSD=$sample-merged-FSD.txt

ls *rmdup.bam > lib.txt

samtools merge $merged `cat lib.txt`

samtools view -f 35 $merged | cut -f 9 | sort | uniq -c | sort -b -k2,2n | sed -e 's/^[/t]*//' > $FSD


