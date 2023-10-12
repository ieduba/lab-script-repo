#! /bin/bash
# SBATCH -N 1
# SBATCH -n 10
echo This is job $SLURM_JOB_ID
echo Test permission

## run this if preseq in the dovetail pipeline doesn't converge. this is the dovetail pip    eline without duplicate removal so picard can run to estimate library complexity.

source activate "microc"

dir=$1 #path to directory with fastqs
sample=$2 #name of fastq up to _R#... (eg K562_low_S72 from K562_low_S72_R1_001.fastq.gz)
fasta=$3 #path to BWA .fasta file (eg /rugpfs/fs0/risc_lab/store/risc_data/downloaded/mm10/genome/Sequence/BWAIndex/mm10_no_alt_analysis_set_ENCODE.fasta)
cores=$4 #number cores to run on

cd $dir
mkdir -p temp

## make input file names
read1=$sample\_R1_001.fastq.gz
read2=$sample\_R2_001.fastq.gz
cut -f1,2 $fasta\.fai > genome.temp

## make output file names
pairs=$sample-wdups-mapped.pairs
bam=$sample-wdups-mapped.PT.bam
picard=$sample-picard.txt

## run two-sided alignment, duplication removal, make pairs file, bam, and stats file
bwa mem -5SP -T0 -t$cores $fasta $read1 $read2 | pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in $cores --nproc-out $cores --chroms-path genome.temp | pairtools sort --tmpdir=temp --nproc $cores | pairtools split --nproc-in $cores --nproc-out $cores --output-pairs $pairs --output-sam -|samtools view -bS -@$cores | samtools sort -@$cores -o $bam
samtools index $bam

## run picard to estimate library complexity
picard EstimateLibraryComplexity I=$bam O=$picard

## housekeeping 
mkdir -p $sample
#mv $sample-* $sample

cd ..
