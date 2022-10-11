#! /bin/bash
# SBATCH -N 1
# SBATCH -n 2 

echo This is job $SLURM_JOB_ID
echo Test permission

source activate "ATACseq"

file=$1
cd LSA2/pipelineRerun

echo 'sample:' $file
cd chromState/$file
pwd

#chromHMM states
for state in {1..15}; do
        echo '  state:' $state
        bamfile=$state-frag_E025_CS.bam
	sortedbam=$state-frag_E025_CS.sort.bam
	fixsortbam=$state-frag_E025_CS.sort.fix.bam
	bedpefile=$state-frag_E025_CS.bedpe
	bedfile=$state-frag_E025_CS.bed
	nrlout=$state-NRL.txt
        samtools sort -n $bamfile > $sortedbam
	samtools fixmate $sortedbam $fixsortbam
	bedtools bamtobed -bedpe -i $fixsortbam > $bedpefile
	awk '{print $1,$2,$6}' $bedpefile > $bedfile
        bash ~/labstore/risc_soft/nrl_finder/nrl_finder -b $bedfile > $nrlout
        bash ~/labstore/risc_soft/nrl_finder/nrl_finder -b $bedfile -o figure 
done

#histone marks and genome wide
for state in H3K9me3 H3K27me3 genomewide; do
        echo '  state:' $state
        bamfile=$state-frag_E025.bam
	sortedbam=$state-frag_E025.sort.bam
	fixsortbam=$state-frag_E025.sort.fix.bam
        bedpefile=$state-frag_E025_CS.bedpe
	bedfile=$state-frag_E025_CS.bed
        nrlout=$state-NRL.txt
	samtools sort -n $bamfile > $sortedbam
	samtools fixmate $sortedbam $fixsortbam
        bedtools bamtobed -bedpe -i $fixsortbam > $bedpefile
	awk '{print $1,$2,$6}' $bedpefile > $bedfile
        bash ~/labstore/risc_soft/nrl_finder/nrl_finder -b $bedfile > $nrlout
        bash ~/labstore/risc_soft/nrl_finder/nrl_finder -b $bedfile -o figure 
done


