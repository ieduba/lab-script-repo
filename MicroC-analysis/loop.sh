#! /bin/bash
# SBATCH -N 1
# SBATCH -n 16

echo This is $SLURM_JOB_ID
echo Test permission

source activate "microc"

hic1=$1 # absolute path to hic file 1
hic2=$2 # absolute path to hic file 2 # SET AS none IF NOT DOING DIFFERENTIAL
celltype=$3 # name of cell type (eg K562 or hela)
diffdir=$4 # path to directory for differential loop results

juicerpath=/rugpfs/fs0/risc_lab/scratch/iduba/linker-histone/Micro-C/juicer
juicertools=$juicerpath/scripts/common/juicer_tools
cores=16

hicdir1=`dirname $hic1`
hicdir2=`dirname $hic2`

echo directory 1: $hicdir1
echo directory 2:$hicdir2

if [ $celltype == K562 ] || [ $celltype == hela ]; then
	genome=hg38
fi
if [ $celltype == mesc ]; then
	genome=mm10
fi

echo cell type: $celltype
echo genome: $genome

for res in 5000 10000 25000; do
	if [[ $hic2 == none ]]; then
		break
	fi
	mkdir -p $diffdir
	cd $diffdir
	mkdir -p $res
	cd ..	
done

for hic in $hic1 $hic2; do
	if [[ $hic == none ]]; then
		break
	fi
	
	hicdir=`dirname $hic`
	cd $hicdir
	pwd
	mkdir -p loops
	cd loops
	
	for res in 5000 10000 25000; do
		mkdir -p $res
		cd $res
		pwd
	
		echo Running loop calling
		bash $juicertools hiccups $hic . -r $res --threads $cores --cpu --ignore_sparsity > hiccups.log 2>&1
		echo Running motif finding
		bash $juicertools motifs $genome $juicerpath/references/motif-$celltype merged_loops.bedpe > motif.log 2>&1
		awk '$25 == "NA" || $25 == "na" || $30 == "NA" || $30 == "na" {print $0}' merged_loops_with_motifs.bedpe > noCTCF_loops.bedpe
		awk '$25 != "NA" && $25 != "na" && $30 != "NA" && $30 != "na" {print $0}' merged_loops_with_motifs.bedpe > CTCF_loops.bedpe
		wc -l *loops.bedpe > Nloop-summary.txt

		echo Running aggregate loop analysis
		for looptype in merged_loops noCTCF_loops CTCF_loops; do
			mkdir -p $looptype-APA
			cd $looptype-APA
### previously ran apa on the same loop set for both hics (both on loops called on hic1), now each runs on its own loop list. which is better?
			bash $juicertools apa -r $res $hic ../$looptype\.bedpe . > apa.log 2>&1
			cp $res/gw/APA.txt ../../../../$diffdir/$res/`basename $hic -mapped-filt-noin-mq30.hic`-$looptype-APA.txt 
			cd ..
		done
		cd ..
	done
	cd ../..
done

for res in 5000 10000 25000; do
	if [[ $hic2 == none ]]; then
		break
	fi

	cd $diffdir
	cd $res
	echo differential loops at resolution $res > diff-loop-summary.txt
	echo Running differential loop analysis at resolution $res
	for looplist in merged_loops.bedpe noCTCF_loops.bedpe CTCF_loops.bedpe; do
		mkdir -p `basename $looplist .bedpe`
		cd `basename $looplist .bedpe`
		bash $juicertools hiccupsdiff $hic1 $hic2 $hicdir1/loops/$res/$looplist $hicdir2/loops/$res/$looplist . --cpu > hiccupsdiff.log 2>&1

		echo $looplist - loops in 1 not in 2 >> ../diff-loop-summary.txt
		wc -l differential_loops1.bedpe >> ../diff-loop-summary.txt
		echo $looplist - loops in 2 not in 1 >> ../diff-loop-summary.txt
		wc -l differential_loops2.bedpe >> ../diff-loop-summary.txt

		echo Plotting differential APA
		for difflooplist in differential_loops1.bedpe differential_loops2.bedpe; do
			diffloopname=`basename $difflooplist .bedpe`-APA
			mkdir -p $diffloopname
			cd $diffloopname
			for hic in $hic1 $hic2; do
				hicname=`basename $hic -mapped-filt-noin-mq30.hic`
				mkdir -p $hicname
				cd $hicname
				bash $juicertools apa -r $res $hic ../../$difflooplist . > apa.log 2>&1
				cp $res/gw/APA.txt ../../$hicname-$diffloopname\.txt
				cd ..
			done
			cd ..
		done
		#file names are hard coded into diffloop python script RIP		
		python ~/scripts/MicroC-analysis/diffloop-apa-normalize.py
		cd ..
	done
	echo Plotting CTCF APA
	python ~/scripts/MicroC-analysis/ctcf-apa-normalize.py
	cd ../..
done
