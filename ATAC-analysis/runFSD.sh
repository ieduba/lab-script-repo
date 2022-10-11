#! /bin/bash
# SBATCH -N 1
# SBATCH -n 2 

echo This is job $SLURM_JOB_ID
echo Test permission

source activate "ATACseq"

file=$1
cd LSA2/pipelineRerun/chromState

echo 'sample:' $file
cd $file
pwd

#chromHMM states
for state in {1..15}; do
        echo '  state:' $state
        CSbam=$state-frag_E025_CS.bam
        samtools index $CSbam
	samtools view -f 35 $CSbam | cut -f 9 | sort | uniq -c | sort -b -k2,2n | sed -e 's/^[ \t]*//' > $state-DIY_FSD.txt
done

#histone marks
echo '  histone marks'
in2=H3K9me3-frag_E025.bam
samtools index $in2
in3=H3K27me3-frag_E025.bam
samtools index $in3
samtools view -f 35 $in2 | cut -f 9 | sort | uniq -c | sort -b -k2,2n | sed -e 's/^[ \t]*//' > H3K9me3-DIY_FSD.txt
samtools view -f 35 $in3 | cut -f 9 | sort | uniq -c | sort -b -k2,2n | sed -e 's/^[ \t]*//' > H3K27me3-DIY_FSD.txt

#genome wide
echo ' genome wide'
in4=genomewide-frag_E025.bam
samtools index $in4
samtools view -f 35 $in4 | cut -f 9 | sort | uniq -c | sort -b -k2,2n | sed -e 's/^[ \t]*//' > genomewide-DIY_FSD.txt



