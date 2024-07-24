#! /bin/bash
#SBATCH -N 1
#SBATCH -n 16

source activate encode-atac-seq-pipeline

infile1=$1 #sorted, filtered bam
infile2=$2 #sorted, filtered bam
compareout=$3

bamCompare -b1 $infile1 -b2 $infile2 -o $compareout -of bigwig --operation ratio -bs 10 -p 16 --scaleFactorsMethod None --normalizeUsing RPKM #scales samples by read count, reports ratio of 1 over 2. bin size 10 bp  
