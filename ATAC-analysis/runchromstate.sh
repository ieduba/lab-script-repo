#! /bin/bash
# SBATCH -N 1
# SBATCH -n 10 

echo This is job $SLURM_JOB_ID
echo Test permission

source activate "nucOcc2"

file=$1

cd ~/labscratch/npagane/pipes/irene_senescence

echo 'sample:' $file
fileshort=`echo $file | cut -d'-' -f1`
fragbam=$fileshort/*.rmdup.bam

#chromHMM states
for state in {1..15}; do
	echo '	state:' $state
	chromstate=~/scripts/LSA2/chromState/CS_bed3s/$state-E025_CS.hg38.bed
	out1=~/scripts/LSA2/pipelineRerun/chromState/$file/$state-frag_E025_CS.bam
	bedtools intersect -a $fragbam -b $chromstate -wa > $out1
done

#histone marks
echo '	histone marks'
out2=~/scripts/LSA2/pipelineRerun/chromState/$file/H3K9me3-frag_E025.bam
out3=~/scripts/LSA2/pipelineRerun/chromState/$file/H3K27me3-frag_E025.bam
bedtools intersect -a $fragbam -b ~/scripts/LSA2/chromState/E025-H3K9me3.hg38.bed -wa > $out2
bedtools intersect -a $fragbam -b ~/scripts/LSA2/chromState/E025-H3K27me3.hg38.bed -wa > $out3

