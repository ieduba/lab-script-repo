#!/bin/bash

#SBATCH --nodes=1
#SBATCH -c 8
#SBATCH -J STAR
#SBATCH -e %x-%j-%a.err
#SBATCH -o %x-%j-%a.out
#SBATCH --verbose

INPUT_R1=$(head -n $SLURM_ARRAY_TASK_ID $inputfiles_R1 | tail -n 1)
INPUT_R2=$(head -n $SLURM_ARRAY_TASK_ID $inputfiles_R2 | tail -n 1)
NAME=$(head -n $SLURM_ARRAY_TASK_ID $samplelist | tail -n 1)
OUTPUT=$outdir\/$NAME\_$SLURM_ARRAY_TASK_ID\_

source activate star-irene

STAR --runThreadN 4 --genomeDir $STAR_INDEX --readFilesIn $INPUT_R1 $INPUT_R2 --outSAMtype BAM SortedByCoordinate \
--outFilterType BySJout --outFilterMultimapNmax 15 --alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 --outFilterMismatchNmax 500 --outFilterMismatchNoverReadLmax 0.05 --alignIntronMin 20 --alignIntronMax 1000000 --outTmpDir $OUTPUT\_tmp \
--outFileNamePrefix $OUTPUT --quantMode GeneCounts --readFilesCommand zcat

