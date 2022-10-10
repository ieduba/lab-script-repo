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

for sample in *mapped-filt.pairs; do
#for sample in */*mESC_E14_12_segments*mapped-filt.pairs; do
	echo $sample
	# separate reads by orientation (inward, outward, or tandem)
	inward=`echo $sample | sed 's/filt/filt-inward/'`
	outward=`echo $sample | sed 's/filt/filt-outward/'`
	tandemp=`echo $sample | sed 's/filt/filt-tandemplus/'`
	tandemm=`echo $sample | sed 's/filt/filt-tandemminus/'`
	awk '$6 == "+" && $7 == "-" { print $0 }' $sample > $inward
	awk '$6 == "-" && $7 == "+" { print $0 }' $sample > $outward
	awk '$6 == "+" && $7 == "+" { print $0 }' $sample > $tandemp
	awk '$6 == "-" && $7 == "-" { print $0 }' $sample > $tandemm

	# adjust fragment start sites so positions are nucleosome centers
#	inadj=$in-adj
#	outadj=$out-adj
#	tandemadj=$sample-tandem-adj
#	awk '{$3+=73; $5-=73}1' $in\.txt > $inadj\.txt
#	awk '{$3-=73; $5+=73}1' $out\.txt > $outadj\.txt
#	awk '{$3+=73; $5+=73}1' $tandemp\.txt > $tandemadj\.txt
#	awk '{$3-=73; $5-=73}1' $tandemm\.txt >> $tandemadj\.txt
	
#	tandem=`echo $sample | sed 's/filt/filt-tandem/'`
#	cat $tandemp > $tandem
#	cat $tandemm >> $tandem

	# subtract positions of cis pairs to get genomic distance
	for file in $inward $outward $tandemp $tandemm; do
		dist=`echo $file | sed 's/.pairs/-dist.txt/'`
		awk '$2 == $4 {print $5 - $3}' $file | sort -n > $dist
	done
done
mv *filt-inward* ../contactProb
mv *filt-outward* ../contactProb
mv *filt-tandem* ../contactProb
#rm *filt-tm*
#rm *filt-tp*
