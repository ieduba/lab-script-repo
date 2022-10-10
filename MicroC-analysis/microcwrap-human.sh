dir=$1 
sample=$2 #fastq name up to _1 or _2

sbatch /rugpfs/fs0/home/iduba/scripts/MicroC-analysis/dovetailQC.sh $dir $sample /rugpfs/fs0/risc_lab/store/risc_data/downloaded/hg38/genome/Sequence/BWAIndex/genome.fa 16 

