#! /bin/bash
# SBATCH -N 1
# SBATCH -n 16

echo This is $SLURM_JOB_ID
echo Test permission

source activate "microc"

hic=$1 # absolute path to hic file 

juicerpath=/rugpfs/fs0/risc_lab/scratch/iduba/linker-histone/Micro-C/juicer
juicertools=$juicerpath/scripts/common/juicer_tools

out=`echo $hic | sed 's/.hic/-domains.txt/'`

## resolution 5000 bp for microc
bash $juicertools arrowhead $hic $out -r 5000 --threads 16 --ignore_sparsity

## resolution 10000 bp (default) for hic
#bash $juicertools arrowhead $hic $out --threads 16 --ignore_sparsity
