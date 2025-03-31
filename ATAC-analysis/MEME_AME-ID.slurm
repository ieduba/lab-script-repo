#! /bin/bash
#SBATCH -N 1
#SBATCH -c 2
#SBATCH 

# example of sbatch command to run this script:
# # sbatch -p risc --array=1-8 --export=backgroundpeakslist=backgroundpeakfiles.txt,inputfiles=inputfiles.txt,outdir=/lustre/fs4/risc_lab/scratch/jyeung/LS_PDvsDoxo_2/ATACseq/post_peakcalling_analysis/04162023/motif_analysis/motif_analysis_output/distalpeaks_in_geneclusters_up,folderlist=folderlist.txt /lustre/fs4/risc_lab/scratch/jyeung/LS_PDvsDoxo_2/ATACseq/post_peakcalling_analysis/04162023/motif_analysis/motif_analysis_input/distalpeaks_in_geneclusters_up/MEME_AME.slurm
# variables:
## outdir=output directory 
## CONTROL=path to control sequences/sequences under your background peak set. 

# this script runs motif enrichment analysis using MEME-AME
INPUT=$(head -n $SLURM_ARRAY_TASK_ID $inputfiles | tail -n 1) # text file containing paths to fasta files
NAME=$(head -n $SLURM_ARRAY_TASK_ID $folderlist | tail -n 1) # text file containing names of future output directories
CONTROL=$(head -n $SLURM_ARRAY_TASK_ID $backgroundpeakslist | tail -n 1) # text file containing paths to fasta files with static peaks (background peak set)
# meme package is installed in Justin's env, JR
source activate JR

# make new directory for each peak set
mkdir $outdir\/$NAME

# run MEME-AME using HOCOMOCOv11_core database
ame --verbose 1 --oc $outdir\/$NAME --scoring avg --method fisher --hit-lo-fraction 0.25 --evalue-report-threshold 10.0 --control $CONTROL --kmer 2 $INPUT /lustre/fs4/risc_lab/store/jyeung/references/motif_databases/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme
