#! /bin/bash
# SBATCH -N 1
# SBATCH -n 4

echo This is job $SLURM_JOB_ID
echo Test permission

source activate rstudio

cd $1
pwd
Rscript pseudorep.R
awk '{print $15,"\t", $2,"\t", $3,"\t",$1,"\t",$6,"\t",$7,"\t",$5}' mnd-1.txt > mnd-1.summary
awk '{print $15,"\t", $2,"\t", $3,"\t",$1,"\t",$6,"\t",$7,"\t",$5}' mnd-2.txt > mnd-2.summary
