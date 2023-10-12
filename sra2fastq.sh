#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16

echo This is job $SLURM_JOB_ID
echo Test permission

source activate data

for SRR in {11142931..11142934}; do
	#prefetch  SRR$SRR
	#fasterq-dump SRR$SRR -t /tmp -e 16 -p
	gzip SRR$SRR/SRR$SRR*.fastq
done

