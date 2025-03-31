#! /bin/bash 
#SBATCH -N 1
#SBATCH -n 4

source activate microc

### run this script to prepare files to use with MEME_AME-ID.slurm. converted from Joanna's R code to bash by Irene

slopdist=50000
genome=~/labstore/risc_data/downloaded/hg38/genome
chroms=$genome/chrom.sizes
fasta=$genome/Sequence/WholeGenomeFasta/genome.fa

mkdir -p MEME-AME-out
rm -f inputfiles.txt controlfiles.txt names.txt

for geneset in up down; do
	# make a bed that's up/down genes +/- 50kb
	bedtools slop -i WTvslow$geneset-genebody.bed -g $chroms -b $slopdist > WTvslow$geneset-genebody-50kb.bed
	for peakset in up down nochange; do
		# find up/down/unchanging ATAC peaks within 50kb of up/down genes
		peakbed=rep-peaks-scrvslow-$peakset\.bed	
		bedtools intersect -a $peakbed -b WTvslow$geneset-genebody-50kb.bed > rep-$peakset-peaks-50kb-$geneset\genes.bed
		# convert beds to fastas for MEME AME motif analysis
		bedtools getfasta -fi $fasta -bed rep-$peakset-peaks-50kb-$geneset\genes.bed > rep-$peakset-peaks-50kb-$geneset\genes.fa
	done

	# write file names to text files that will be inputs for MEME script. unchanigng peaks are control for both up and down peaks.
	echo rep-up-peaks-50kb-$geneset\genes.fa >> inputfiles.txt
	echo rep-down-peaks-50kb-$geneset\genes.fa >> inputfiles.txt
	echo rep-nochange-peaks-50kb-$geneset\genes.fa >> controlfiles.txt
	echo rep-nochange-peaks-50kb-$geneset\genes.fa >> controlfiles.txt
	echo up-peaks-$geneset-genes >> names.txt
	echo down-peaks-$geneset-genes >> names.txt
done
