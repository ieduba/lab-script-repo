#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16

echo This is job $SLURM_JOB_ID
echo Test permission

source activate data

for SRR in {10560118..10560129}; do
#	prefetch  SRR$SRR
#	fasterq-dump SRR$SRR -t /tmp -e 16 -p
	gzip SRR$SRR*.fastq
done

