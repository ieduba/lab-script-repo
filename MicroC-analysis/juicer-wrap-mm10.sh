#! /bin/bash
# SBATCH -N 1
# SBATCH -n 8

echo This is $SLURM_JOB_ID
echo Test permission

source activate "microc"

bash /rugpfs/fs0/risc_lab/scratch/iduba/linker-histone/Micro-C/juicer/scripts/juicer.sh -d . -z /rugpfs/fs0/risc_lab/store/risc_data/downloaded/mm10/genome/Sequence/BWAIndex/mm10_no_alt_analysis_set_ENCODE.fasta -p /rugpfs/fs0/risc_lab/store/risc_data/downloaded/mm10/genome/chrom.sizes -a $1 -D /rugpfs/fs0/risc_lab/scratch/iduba/linker-histone/Micro-C/juicer 
