#! /bin/bash
# SBATCH -N 1
# SBATCH -n 4

source activate encode-atac-seq-pipeline

# try to match encode peak calling as closely as possible: https://github.com/ENCODE-DCC/chip-seq-pipeline2/releases/tag/v1.3.5.1 but with gapped peaks (the result of calling broad peaks. min-length and max-gap can be used to tune gapped peak behavior but not sure what makes sense)

chipbam=$1
inbam=$2
prefix=`echo $chipbam | sed 's/.bam//'`
outdir='.'

macs2 callpeak --broad -t $chipbam -c $inbam -f BAM -n $prefix -g hs --outdir $outdir --keep-dup all -p 0.01
