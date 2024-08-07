#! /bin/bash
# SBATCH -N 1
# SBATCH -n 30
echo This is job $SLURM_JOB_ID
echo Test permission

## dovetail_qc.sh by Irene Duba October 2021. Runs dovetail recommended QC on micro-c data (zipped fastqs): makes .pairs and .bam files, tabulates length 
## distributions, runs preseq to estimate library size

source activate "microc"

dir=$1 #path to directory with fastqs
sample=$2 #name of fastq up to _R#... (eg K562_low_S72 from K562_low_S72_R1_001.fastq.gz)
fasta=$3 #path to BWA .fasta file (eg /rugpfs/fs0/risc_lab/store/risc_data/downloaded/mm10/genome/Sequence/BWAIndex/mm10_no_alt_analysis_set_ENCODE.fasta)
cores=$4 #number cores to run on

cd $dir
mkdir -p temp

## make input file names
read1=$sample\_R1.fastq.gz
read2=$sample\_R2.fastq.gz
cut -f1,2 $fasta\.fai > genome.temp

## make output file names
stats=$sample-stats.txt
pairs=$sample-mapped.pairs
bam=$sample-mapped.PT.bam
qc=$sample-QC.txt
pseq=$sample-preseq.txt

## run two-sided alignment, duplication removal, make pairs file, bam, and stats file. separate steps to be able to pick up mid-way if it fails
bwa mem -5SP -T0 -t$cores $fasta $read1 $read2 -o $sample-aligned.sam
pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in $cores --nproc-out $cores --chroms-path genome.temp $sample-aligned.sam > $sample-parsed.pairsam
pairtools sort --tmpdir=temp --nproc $cores $sample-parsed.pairsam > $sample-sorted.pairsam
pairtools dedup --nproc-in $cores --nproc-out $cores --mark-dups --output-stats $stats --output $sample-dedup.pairsam $sample-sorted.pairsam
pairtools split --nproc-in $cores --nproc-out $cores --output-pairs $pairs --output-sam $sample-unsorted.bam $sample-dedup.pairsam
samtools sort -@$cores -o $bam $sample-unsorted.bam
samtools index $bam

## run python QC script to simplify stats file
python /rugpfs/fs0/risc_lab/scratch/iduba/linker-histone/Micro-C/get_qc_edit.py -p $stats > $qc 

## housekeeping 
#rm genome.temp
# rm -r temp
mkdir -p $sample
mv $sample-* $sample
mkdir $sample/intermediates
mv $sample-aligned.sam $sample-parsed.pairsam $sample-sorted.pairsam $sample-dedup.pairsam $sample-unsorted.bam $sample/intermediates


cd ..
