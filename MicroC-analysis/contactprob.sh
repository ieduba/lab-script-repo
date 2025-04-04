#! /bin/bash
# SBATCH -N 1
# SBATCH -n 8

echo This is job $SLURM_JOB_ID
echo Test permission

## contactprob.sh by Irene Duba August 2021.
## calculates genomic distance between Micro-C pairs, adjusted for ligation orientation

mkdir -p contactProb
cd $1
pwd

#for pair in *-mapped-filt.pairs; do #may want to change what text it searches for to run on subset of pairs files
for pair in *allTss-ENCFF963KIA-mapped-filt.pairs; do
	echo $pair
	sample=$pair	

## uncomment below if not yet filtered (anno-intersect.sh script does this filtering step)
	#sample=`echo $pair | sed 's/mapped/mapped-filt/'`
	#awk '$2 !~ /chrM|chrY|chrU.+|chr.+_.+/ && $4 !~ /chrM|chrY|chrU.+|chr.+_.+/ {print $0}' $pair > $sample

	# separate reads by orientation (inward, outward, or tandem)
	inward=`echo $sample | sed 's/mapped-filt/mapped-filt-inward/'`
	outward=`echo $sample | sed 's/mapped-filt/mapped-filt-outward/'`
	tandemp=`echo $sample | sed 's/mapped-filt/mapped-filt-tandemplus/'`
	tandemm=`echo $sample | sed 's/mapped-filt/mapped-filt-tandemminus/'`
	awk '$6 == "+" && $7 == "-" { print $0 }' $sample > $inward
	awk '$6 == "-" && $7 == "+" { print $0 }' $sample > $outward
	awk '$6 == "+" && $7 == "+" { print $0 }' $sample > $tandemp
	awk '$6 == "-" && $7 == "-" { print $0 }' $sample > $tandemm

	# subtract positions of cis pairs (start of each read) to get genomic distance
	for file in $inward $outward $tandemp $tandemm; do
		dist=`echo $file | sed 's/.pairs/-dist.txt/'`
		awk '$2 == $4 {print $5 - $3}' $file | sort -n > $dist
	done
done
mv *filt-inward*txt ../contactProb
mv *filt-outward*txt ../contactProb
mv *filt-tandem*txt ../contactProb
