#! /bin/bash
# SBATCH -N 1
# SBATCH -n 32

echo This is $SLURM_JOB_ID
echo Test permission

source activate "microc"

bash /rugpfs/fs0/risc_lab/scratch/iduba/linker-histone/Micro-C/juicer/scripts/juicer.sh -d . -z /rugpfs/fs0/risc_lab/store/risc_data/downloaded/hg38/genome/Sequence/BWAIndex/genome.fa -p /rugpfs/fs0/risc_lab/store/risc_data/downloaded/hg38/genome/chrom.sizes -a $1 -D /rugpfs/fs0/risc_lab/scratch/iduba/linker-histone/Micro-C/juicer -t 32 #-S chimeric 

#bash /rugpfs/fs0/risc_lab/scratch/iduba/linker-histone/Micro-C/juicer/scripts/juicer.sh -d . -g hg38 -a $1 -D /rugpfs/fs0/risc_lab/scratch/iduba/linker-histone/Mirco-C/juicer
