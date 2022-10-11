#! /bin/bash
# SBATCH -N 1
# SBATCH -N 1

echo This is job $SLURM_JOB_ID
echo Test permission

#source activate "ATACseq"

file=$1
sorted=`basename $file | sed ' s/.trim.st.all.blft.qft.bam/-sorted.bam/'`	
#samtools sort $file > $sorted 

#conda deactivate
#source activate "senescence"

out=`basename $sorted | sed ' s/sorted.bam/extrap.txt/'`
#preseq lc_extrap -B -P -o $out $sorted

##find max in second column of preseq output (max # unique reads) and divide by 2 to set the threshold. threshold for total reads is in the first column.
thresh=`awk -v m=0 '{if ($2>0+m) m=$2; t=(m/2)} END{print t}' $out`
threshTot=`awk -v tU="$thresh" '{if($2<tU){threshTot=$1}}END{print threshTot}' $out `

echo -e $out "\t" $threshTot "\t" $thresh >> thresholds.txt
