#! /bin/bash
#SBATCH -N 1
#SBATCH -n 16

source activate encode-atac-seq-pipeline

infile1=$1 #sorted, filtered bam
outfile1=`echo $infile1 | sed 's/.bam/.bed/'`

infile2=$2 #sorted, filtered bam
outfile2=`echo $infile2 | sed 's/.bam/.bed/'`

compareout=$3

bamCoverage -b $infile1 -o $outfile1 --outFileFormat bedgraph #bin size is 50 bp
bamCoverage -b $infile2 -o $outfile2 --outFileFormat bedgraph #bin size is 50 bp

bamCompare -b1 $infile1 -b2 $infile2 -o $compareout --outFileFormat bedgraph --operation ratio #uses read count to normalize between samples, reports ratio of 1 over 2. bin size 50 bp  
