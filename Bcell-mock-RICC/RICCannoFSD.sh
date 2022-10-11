#! /bin/bash
# SBATCH -N 1
# SBATCH -n 2 

echo This is job $SLURM_JOB_ID
echo Test permission

source activate "ATACseq"

directory=$1
cd ~/linker-histone/Bcell/
cd $directory
pwd

echo 'controls'
controls=`ls *controlbins*bam`
for ctrl in $controls; do
	echo $ctrl
	shortname=`echo $ctrl | cut -d'.' -f1`
	out=$ctrl-FSD.txt
	samtools index $ctrl
	samtools view -f 35 $ctrl | cut -f 9 | sort | uniq -c | sort -b -k2,2n | sed -e 's/^[/t]*//' > $out
done

echo 'ROIs'
ROIs=`ls hg19*BJ*bam`
for ROI in $ROIs; do
	echo $ROI
	shortname=`echo $ctrl | cut -d'.' -f1`
	out=$ROI-FSD.txt
	samtools index $ROI
	samtools view -f 35 $ROI | cut -f 9 | sort | uniq -c | sort -b -k2,2n | sed -e 's/^[/t]*//' > $out
done

