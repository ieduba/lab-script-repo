#!/bin/bash

#SBATCH --nodes=1
#SBATCH -c 8
#SBATCH -J UMI
#SBATCH -e %x-%j-%a.err
#SBATCH -o %x-%j-%a.out
#SBATCH --verbose

INPUT_R1=$(head -n $SLURM_ARRAY_TASK_ID $inputfiles_R1 | tail -n 1)
INPUT_R2=$(head -n $SLURM_ARRAY_TASK_ID $inputfiles_R2 | tail -n 1)
NAME=$(head -n $SLURM_ARRAY_TASK_ID $samplelist | tail -n 1)
OUTPUT=$outdir\/$NAME\_$SLURM_ARRAY_TASK_ID\_

source activate star-irene

fastp -i $INPUT_R1 -I $INPUT_R2 -o $OUTPUT-umi_R1.fastq.gz -O $OUTPUT-umi_R2.fastq.gz -U --umi_loc per_read --umi_skip 2 --umi_len 5 -D -h $OUTPUT\.html 

