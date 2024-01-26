#! /bin/bash
# SBATCH -N 1
# SBATCH -n 4

source activate encode-atac-seq-pipeline

chipbam='scr-3-H1/scr-3-H1_S51_001.trim.st.all.blft.qft.rmdup.bam'
inputbam='scr-3-input/scr-3-input_S50_001.trim.st.all.blft.qft.rmdup.bam'
prefix='scr-H1-pilot-deeper'
outdir='broadpeaks'

macs2 callpeak --broad -t $chipbam -c $inputbam -f BAMPE -n $prefix -g hs --outdir $outdir 
