#!/bin/bash
#SBATCH -N 1
#SBATCH -n 10

echo This is job $SLURM_JOB_ID
echo Test permission

#fetch SRA files from online
for SRR in {8954797..8954830}; do
	prefetch  SRR$SRR
done

#convert .sra files to fastqs
for sra in `ls SRR*/*.sra`; do
	fasterq-dump $sra
done

#zip fastqs 
for fastq in `ls *.fastq`; do
	gzip $fastq
done 
