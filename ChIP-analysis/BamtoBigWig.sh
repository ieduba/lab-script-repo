#! /bin/bash
#SBATCH -N 1
#SBATCH -n 16

source activate encode-atac-seq-pipeline

infile1=$1 #sorted, filtered bam
outfile1=`echo $infile1 | sed 's/.bam/.bed/'`

bamCoverage -b $infile1 -o $outfile1 -bs 10 -of bedgraph --normalizeUsing RPGC --effectiveGenomeSize 2913022398 #normalized by reads per genome coverage, effective genome size from deeptoolsdocumentation

## EFFECTIVE GENOME SIZE: hg38: 2913022398, mm10: 2652783500

