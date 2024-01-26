#! /bin/bash
#SBATCH -N 1
#SBATCH -n 2

source activate microc

mkdir -p intersections
region_dir='/ru-auth/local/home/iduba/linker-histone/Micro-C/annotations/chip/K562'
for bed in *.bed; do
	for region in H3K9me3 H3K27ac H3K27me3; do
		region_path=$region_dir/K562-$region\.st.bed
		out=intersections/$region-$bed
		bedtools intersect -a $bed -b $region_path -wa > $out
	done
done 
