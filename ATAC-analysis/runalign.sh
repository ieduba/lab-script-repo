#! /bin/bash
# SBATCH -N 1
# SBATCH -n 20 

echo This is job $SLURM_JOB_ID
echo Test permission

source activate "ATACseq-senescence"

srun hg19align.sh . libs.txt 


