#! /bin/bash
# SBATCH -N 1
# SBATCH -n 4

echo This is job $SLURM_JOB_ID
echo Test permission

source activate microc

sample=$1 #name of sample 
bedin=$2 #name of bed file to split into chunks
Nchunks=$3 #up to 100

echo $sample

bam=$sample-mapped.PT.bam
pair=$sample-mapped.pairs

#totallines=`wc -l < $bedin`
#chunksize=$(($totallines/$Nchunks))
name=`basename $bedin | sed 's/.bed/_/'`

#split -d -l $chunksize --verbose $bedin $name 

for chunk in `ls $name??`; do
	echo making bams for $chunk
	chunkbam=$chunk-$bam
	chunkpair=$chunk-$pair
	# make new bams with only reads that fall in the chunk 
	bedtools intersect -a $bam -b $chunk -wa > $chunkbam
	samtools view $chunkbam > $chunk-temp.sam

	# from euvshet to remove lines with names that show up in both bams - not sure how to implement for chunks
	#awk 'NR==FNR{A[$1]=$2;next} !A[$1]{print}' temphet.sam tempeu.sam > tempeu2.sam
	#awk 'NR==FNR{A[$1]=$2;next} !A[$1]{print}' tempeu.sam temphet.sam > temphet2.sam

	echo making pairs for $chunk
	# filter .pairs files to include only lines from  eu or het bams (matched by names) to make
	#       eu and het .pairs files for contact probability analysis
	awk 'NR==FNR{A[$1]=$2;next} A[$1]{print}' $chunk-temp.sam $pair > $chunkpair
done
rm *temp.sam

