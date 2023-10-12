dir=$1 

for sample in hesc_BR1_TR1_1 hesc_BR1_TR1_2 hesc_BR1_TR2_1 hesc_BR1_TR2_2; do
	sbatch -p risc,hpc /rugpfs/fs0/home/iduba/scripts/MicroC-analysis/dovetailQC.sh $dir $sample /rugpfs/fs0/risc_lab/store/risc_data/downloaded/hg38/genome/Sequence/BWAIndex/genome.fa 10 
done
