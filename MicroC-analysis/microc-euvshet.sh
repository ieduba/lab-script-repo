#! /bin/bash
# SBATCH -N 1
# SBATCH -n 8

echo This is job $SLURM_JOB_ID
echo Test permission

source activate microc

# make file names. zippedbedin is ChromHMM bed file with 12 states E1 - E12 or 18 named states..
bedname=$1
#echo $bedin
sample=$2
echo $sample
#bedname=`basename $bedin`
echo $bedname
eubed=eu-$bedname
hetbed=het-$bedname
#state=$3

#echo splitting bed
#if [ $state -eq 12 ]; then
#	# assign lines of the ChromHMM bed file to heterochromatin or euchromatin bed files as follows:
#	# Heterochromatin: 1 (CTCF "insulator"), 2 (intergenic region), 3 (H3K9me3 "heterochromatin"), 
#	#	5 (H3K27me3 "repressed chromatin"), 12 (low signal/repeat)
#	# Euchromatin: 4 (enhancer), 7 (active promoter), 8 (strong enhancer), 9 (transcriptional transition), 
#	#	10 (transcriptional elongation), 11 (weak/poised enhancer) 
#	# 6 (bivalent promoter) is not assigned to either.
##	gunzip $zippedbedin
#	awk '$4 == "E1" || $4 == "E2" || $4 == "E3" || $4 == "E5" || $4 == "E12" {print $1 "\t" $2 "\t" $3}' $bedin > $hetbed
#	awk '$4 == "E4" || $4 == "E7" || $4 == "E8" || $4 == "E9" || $4 == "E10" || $4 == "E11" {print $1 "\t" $2 "\t" $3}' $bedin > $eubed
##	gzip $bedin
#fi
#
#if [ $state -eq 18 ]; then
	# assign lines of the ChromHMM  bed file to heterochromatin or euchromatin bed files using names
#	awk '$4 == "TssA" || $4 == "TssFlnk" || $4 == "TssFlnkU" || $4 == "TssFlnkD" || $4 == "Tx"|| $4 == "TxWk" || $4 == "EnhG1" || $4 == "EnhG2" || $4 == "EnhA1" || $4 == "EnhA2" || $4 == "EnhWk" {print $1 "\t" $2 "\t" $3}' $bedin > $eubed
#	awk '$4 == "Het" || $4 == "ReprPC" || $4 == "ReprPCWk" || $4 == "Quies" {print $1 "\t" $2 "\t" $3}' $bedin > $hetbed
#fi

# making more file names
bam=$sample-mapped.PT.bam
hetbam=het-$bam
eubam=eu-$bam
pair=$sample-mapped-filt.pairs

echo intersecting bams
# make new bams with only reads that fall in heterochromatin or euchromatin regions 
bedtools intersect -a $bam -b $hetbed -wa > $hetbam
bedtools intersect -a $bam -b $eubed -wa > $eubam

echo making sams
# turn bams to sams for awk manipluation
samtools view $hetbam > $sample-temphet.sam
samtools view $eubam > $sample-tempeu.sam

# remove lines with names that show up in both eu and het bams
awk 'NR==FNR{A[$1]=$2;next} !A[$1]{print}' $sample-temphet.sam $sample-tempeu.sam > $sample-tempeu2.sam
awk 'NR==FNR{A[$1]=$2;next} !A[$1]{print}' $sample-tempeu.sam $sample-temphet.sam > $sample-temphet2.sam

echo filtering pairs
# filter .pairs files to include only lines from  eu or het bams (matched by names) to make
# 	eu and het .pairs files for contact probability analysis
awk 'NR==FNR{A[$1]=$2;next} A[$1]{print}' $sample-temphet2.sam $pair > het-$pair
awk 'NR==FNR{A[$1]=$2;next} A[$1]{print}' $sample-tempeu2.sam $pair > eu-$pair

rm $sample-temphet.sam $sample-tempeu.sam $sample-temphet2.sam $sample-tempeu2.sam

