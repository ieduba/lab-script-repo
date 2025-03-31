#! /bin/bash
#SBATCH -N 1
#SBATCH -n 4

juicerpath=/rugpfs/fs0/risc_lab/scratch/iduba/linker-histone/Micro-C/juicer
juicertools=$juicerpath/scripts/common/juicer_tools

for TSS in WTvslowup-loops.bed WTvslow-nochange-loops.bed; do
	TSSname=`basename $TSS -loops.bed`-APA
	mkdir -p $TSSname
	cd $TSSname
	for hic in scr-combo-mapped-filt-noin-mq30.hic dH1-combo-mapped-filt-noin-mq30.hic; do
		hicname=`basename $hic -mapped-filt-noin-mq30.hic`
		mkdir -p $hicname
		cd $hicname
		bash $juicertools apa -r 5000 ../../../$hic ../../$TSS . > apa.log 2>&1
		cp 5000/gw/APA.txt ../../$hicname-$TSSname\.txt
		cd ..
	done
	cd ..
done
#file names are hard coded into diffgene python script RIP              
python diffgene-apa-normalize.py

for TSS in K562_scr1_q1-loops.bed K562_scr1_q4sample-loops.bed; do
	TSSname=`basename $TSS -loops.bed`-APA
	mkdir -p $TSSname
	cd $TSSname
	for hic in scr-combo-mapped-filt-noin-mq30.hic dH1-combo-mapped-filt-noin-mq30.hic; do
		hicname=`basename $hic -mapped-filt-noin-mq30.hic`
		mkdir -p $hicname
		cd $hicname
		bash $juicertools apa -r 5000 ../../../$hic ../../$TSS . > apa.log 2>&1
		cp 5000/gw/APA.txt ../../$hicname-$TSSname\.txt
		cd ..
	done
	cd ..
done
#file names are hard coded into geneqs python script RIP              
python geneqs-apa-normalize.py
