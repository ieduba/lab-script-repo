#! /bin/bash
# SBATCH -N 1
# SBATCH -n 24

echo This is job $SLURM_JOB_ID
echo Test permission

res=$1 ## resolution in bp
sample=$2

cd $sample/$sample
pwd
for chrom in {1..22}; do
	java -jar  /rugpfs/fs0/risc_lab/scratch/iduba/linker-histone/Micro-C/juicer/CPU/common/juicer_tools.1.9.9_jcuda.0.8.jar eigenvector -p KR $sample-mapped-filt-mq30.hic $chrom BP $res $sample\_$chrom\_eigen_$res\.txt
done
cd ../..
