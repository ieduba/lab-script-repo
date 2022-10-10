#! /bin/bash
# SBATCH -N 1
# SBATCH -n 8

echo This is job $SLURM_JOB_ID
echo Test permission

for sample in mESC-WT mESC-low; do
	cd $sample/aligned
	pwd
	for chrom in {1..19}; do
		java -jar  /rugpfs/fs0/risc_lab/scratch/iduba/linker-histone/Micro-C/juicer/CPU/common/juicer_tools.1.9.9_jcuda.0.8.jar eigenvector -p KR inter_30.hic $chrom BP 100000 $sample\_$chrom\_eigen.txt
	done
	cd ../..
done
