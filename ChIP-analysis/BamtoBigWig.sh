#! /bin/bash
#SBATCH -N 1
#SBATCH -n 16

source activate encode-atac-seq-pipeline

infile1=$1 #sorted, filtered bam
outfile1=`echo $infile1 | sed 's/.bam/.bw/'`

bamCoverage -b $infile1 -o $outfile1 -bs 10 --normalizeUsing RPKM  #normalized by reads per kilobase per million mapped reads 

