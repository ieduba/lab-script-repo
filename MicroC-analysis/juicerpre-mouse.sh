#! /bin/bash
# SBATCH -N 1
# SBATCH -n 8

echo This is job $SLURM_JOB_ID
echo Test permission

for sample in mESC-WT mESC-low; do
	cd $sample/aligned
	pwd
	for rep in saddle-1.txt 1 2; do
		java -jar  /rugpfs/fs0/risc_lab/scratch/iduba/linker-histone/Micro-C/juicer/CPU/common/juicer_tools.1.9.9_jcuda.0.8.jar pre saddle-$rep\.txt saddle-$rep\.hic /rugpfs/fs0/risc_lab/store/risc_data/downloaded/mm10/genome/chrom.sizes
	done
	cd ../..
done
