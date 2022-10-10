#! /bin/bash
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --verbose
#SBATCH -J pseudo-RNA
#SBATCH -e %x-%j-%a.err
#SBATCH -o %x-%j-%a.out

INPUT_R1=$(head -n $SLURM_ARRAY_TASK_ID $inputfiles_R1 | tail -n 1)
INPUT_R2=$(head -n $SLURM_ARRAY_TASK_ID $inputfiles_R2 | tail -n 1)
NAME=$(head -n $SLURM_ARRAY_TASK_ID $samplelist | tail -n 1)
OUTPUT=$outdir\/$NAME

cd $outdir

conda activate jr_rna

seqtk sample -s2022 $INPUT_R1 0.5 > $OUTPUT\_1_pr1.fastq.gz
seqtk sample -s2022 $INPUT_R2 0.5 > $OUTPUT\_2_pr1.fastq.gz
seqtk sample -s23 $INPUT_R1 0.5 > $OUTPUT\_1_pr2.fastq.gz
seqtk sample -s23 $INPUT_R2 0.5 > $OUTPUT\_2_pr2.fastq.gz

kallisto index -i transcripts.idx $transcript.fasta

kallisto quant -i transcripts.idx -o $OUTPUT\_pr1_kallisto -b 100 <(zcat $OUTPUT\_1_pr1.fastq.gz) <(zcat $OUTPUT\_2_pr1.fastq.gz)
kallisto quant -i transcripts.idx -o $OUTPUT\_pr2_kallisto -b 100 <(zcat $OUTPUT\_1_pr2.fastq.gz) <(zcat $OUTPUT\_2_pr2.fastq.gz)
