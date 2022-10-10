dir=$1
sample=$2

sbatch /rugpfs/fs0/home/iduba/scripts/MicroC-analysis/dovetailQC.sh $dir $sample /rugpfs/fs0/risc_lab/store/risc_data/downloaded/mm10/genome/Sequence/BWAIndex/mm10_no_alt_analysis_set_ENCODE.fasta 16

