#! /bin/bash
#SBATCH -N 1
#SBATCH -n 1

source activate microc

bed1=`echo $1 | sed 's/.bedpe/.bed/'`
bed2=`echo $2 | sed 's/.bedpe/.bed/'`

cut -f 1-3,12 $1 | tail -n +3 > $bed1
cut -f 1-3,12 $2 | tail -n +3 > $bed2

namepart1=`echo $bed1 | sed 's/.bed//'`
name=$namepart1-$bed2

bedtools intersect -wo -f 0.9 -r -a $bed1 -b $bed2 > $name
