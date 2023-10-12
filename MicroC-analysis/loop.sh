#! /bin/bash
# SBATCH -N 1
# SBATCH -n 16

echo This is $SLURM_JOB_ID
echo Test permission

source activate "microc"

hic1=$1 # absolute path to hic file 1
hic2=$2 # absolute path to hic file 2 # SET AS none IF NOT DOING DIFFERENTIAL
celltype=$3 # name of cell type (eg K562 or hela)
diffdir=$4 # absolute path to directory for differential loop results
#loops=/rugpfs/fs0/risc_lab/scratch/iduba/linker-histone/Micro-C/annotations/$genome-TSS-TTS.bedpe

juicerpath=/rugpfs/fs0/risc_lab/scratch/iduba/linker-histone/Micro-C/juicer
juicertools=$juicerpath/scripts/common/juicer_tools
cores=16

hicdir1=`dirname $hic1`
hicdir2=`dirname $hic2`

#if [[$celltype == K562]] || [[$celltype == hela]]; then
#	genome=hg38
#fi
#if [[$celltype == mesc]]; then
#	genome=mm10
#fi

genome=mm10
echo $celltype
echo $genome

for hic in $hic1 $hic2; do
	if [[ $hic == none ]]; then
		break
	fi
	
	hicdir=`dirname $hic`
	cd $hicdir
	pwd
#	mkdir -p TSS-TTS
#	cd TSS-TTS
	mkdir -p loops
	cd loops
	
	for res in 5000 10000 25000; do
		mkdir -p $res
		cd $res
		pwd
	
		echo Running loop calling
# add $loops right before -r for loop list
		bash $juicertools hiccups $hic . -r $res --threads $cores --cpu --ignore_sparsity > hiccups.log 2>&1
		echo Running motif finding
		bash $juicertools motifs $genome $juicerpath/references/motif-$celltype merged_loops.bedpe > motif.log 2>&1
		awk '$25 == "NA" || $25 == "na" {print $0}' merged_loops_with_motifs.bedpe > noCTCF-loops.bedpe
		awk '$25 != "NA" && $25 != "na" {print $0}' merged_loops_with_motifs.bedpe > CTCF-loops.bedpe
		wc -l *loops.bedpe > Nloop-summary.txt

# runs on both using the motif annotated loop lists from the first input
		echo Running aggregate loop analysis
#		for looplist in merged_loops.bedpe; do
		for looplist in merged_loops.bedpe noCTCF-loops.bedpe CTCF-loops.bedpe; do
			mkdir -p `basename $looplist .bedpe`-APA
			cd `basename $looplist .bedpe`-APA
			bash $juicertools apa -r 10000 $hic $hicdir1/loops/$res/$looplist . > apa.log 2>&1
			cd ..
		done
		cd ..
	done
done

for res in 5000 10000 25000; do
	if [[ $hic2 == none ]]; then
		break
	fi

	cd $diffdir
	mkdir -p $res
	cd $res
	echo differential loops at resolution $res > diff-loop-summary.txt

	echo Running differential loop analysis at resolution $res
	for looplist in merged_loops.bedpe noCTCF-loops.bedpe CTCF-loops.bedpe; do
#	for looplist in merged_loops.bedpe; do
		mkdir -p `basename $looplist .bedpe`
		cd `basename $looplist .bedpe`
		bash $juicertools hiccupsdiff $hic1 $hic2 $hicdir1/loops/$res/$looplist $hicdir2/loops/$res/$looplist . --cpu > hiccupsdiff.log 2>&1
		echo $looplist - loops in 1 not in 2 >> ../diff-loop-summary.txt
		wc -l differential_loops1.bedpe >> ../diff-loop-summary.txt
		echo $looplist - loops in 2 not in 1 >> ../diff-loop-summary.txt
		wc -l differential_loops2.bedpe >> ../diff-loop-summary.txt
		cd ..
	done
done




