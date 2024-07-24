#! /bin/bash
#SBATCH -N 1
#SBATCH -n 4

cspath=~/scripts/CscoreTool

pairs=$1
echo $pairs
summary=`echo $pairs | sed 's/.pairs/-cis-5kb.summary/'`
prefix=`echo $pairs | sed 's/.pairs/-cis-5kb/'`

echo making $summary
#double check how much you're tailing -- seems like the length of the pairs header varies?
tail -n +401 $pairs | awk -v OFS="\t" '($2==$4)&&($5-$3+1>5000){print $1,$2,$3,$6,$4,$5,$7}' > $summary

echo  running $cspath/CscoreTool1.1 $cspath/hg38-100kb.bed $summary $prefix 4 1000000
$cspath/CscoreTool1.1 $cspath/hg38-100kb.bed $summary $prefix 4 1000000
