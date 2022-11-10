#! /bin/bash
# SBATCH -N 1
# SBATCH -n 6

echo This is job $SLURM_JOB_ID
echo Test permission

cd $1
pwd

for pair in */*mapped.pairs; do
	echo $pair
	filtered=`echo $pair | sed 's/mapped/mapped-filt/'`
	awk '$2 !~ /chrM|chrY|chrU.+|chr.+_.+/ && $4 !~ /chrM|chrY|chrU.+|chr.+_.+/ {print $0}' $pair > $filtered
done
