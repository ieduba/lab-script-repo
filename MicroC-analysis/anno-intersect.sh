#! /bin/bash
#SBATCH -N 1
#SBATCH -n 2

source activate microc

sample=$1
bed=$2
echo $bed
bam=$sample-mapped.PT.bam
echo $bam
rawpair=$sample-mapped.pairs
echo $rawpair
anno=`basename $bed .bed`

# filter pairs
pair=`echo $rawpair | sed 's/mapped/mapped-filt/'`
awk '$2 !~ /chrM|chrY|chrU.+|chr.+_.+/ && $4 !~ /chrM|chrY|chrU.+|chr.+_.+/ {print $0}' $rawpair > $pair

# make new bams with only reads that fall in region 
bedtools intersect -a $bam -b $bed -wa | samtools view > $anno-temp.sam

# filter .pairs files to include only lines from anno bam (matched by names) to make
#       anno .pairs file for contact probability analysis
awk 'NR==FNR{A[$1]=$2;next} A[$1]{print}' $anno-temp.sam $pair > $sample-$anno-mapped-filt.pairs

