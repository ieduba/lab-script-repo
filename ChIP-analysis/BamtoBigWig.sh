#! /bin/bash
#SBATCH -N 1
#SBATCH -n 4

source activate encode-atac-seq-pipeline

infile1=$1 #filtered bam
sorted=`echo $infile1 | sed 's/.bam/.st.bam/'`
outfile1=`echo $infile1 | sed 's/.bam/.bw/'`

samtools sort $infile1 > $sorted
samtools index $sorted
bamCoverage -b $sorted -o $outfile1 -bs 1000 -of bigwig --normalizeUsing RPGC --effectiveGenomeSize 2913022398 #normalized by reads per genome coverage, effective genome size from deeptools documentation

## EFFECTIVE GENOME SIZE: hg38: 2913022398, mm10: 2652783500

