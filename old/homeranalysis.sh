#! /bin/bash
# SBATCH -N 1
# SBATCH -n 10 

echo This is job $SLURM_JOB_ID
echo Test permission

source activate "homer-analysis"

beddir=$1
genome=$2

if [ -z "$genome" ]
then
	genome=/rugpfs/fs0/risc_lab/store/risc_data/downloaded/hg19/genome/Sequence/Bowtie2Index/genome.fa
fi

cd $beddir

for comp in CyclingvsQuiescent CyclingvsSen_09d CyclingvsSen_21d QuiescentvsSen_09d QuiescentvsSen_21d Sen_09dvsSen_21d; do
	for level in Down Up; do
		bedin=$level${comp}.bed	
		mkdir HomerOut_$level$comp
		outdir=HomerOut_$level$comp
		background=Norm${comp}.bed
		findMotifsGenome.pl $bedin $genome $outdir -size 200 -bg $background -nomotif
	done
	# find known motifs in random sampling of peaks (for p-value calibration)
	randbed=${comp}_randFG.bed
	randBGbed=${comp}_randBG.bed
	mkdir HomerOut_rand$comp
	outdir=HomerOut_rand$comp
	findMotifsGenome.pl $randbed $genome $outdir -size 200 -bg $randBGbed -nomotif
done

