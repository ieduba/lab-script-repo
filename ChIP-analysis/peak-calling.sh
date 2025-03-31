#! /bin/bash
#SBATCH -N 1
#SBATCH -n 8

source activate cut_tag
 
input=$1 #filtered bam from fastq2bam pipeline
control=$2 #IgG or other control, filtered bam from fastq2bam pipeline

out_dir=`basename $input _001.trim.st.all.blft.qft.rmdup.bam`
mkdir -p $out_dir

#gap=2000 #Arnold found that anything bigger than this fails, we use bedtools below to make 10kb gapped peaks
gap=600

sicer -t $input -c $control -s hg38 -o . -fdr 0.1 -g $gap -cpu 8
mv $out_dir*island* $out_dir

##merge peaks within 10kb of each other to artificially get gapped peaks with 10kb gap
#bed=$out_dir/*island.bed
#stbed=`echo $bed | sed 's/island.bed/island.st.bed/'`
#cat $bed | sed 's/\t$//' | bedtools sort -i - > $stbed #make bed tab delimited before sorting
#bedtools merge -i $stbed -d 10000 > $out_dir/$out_dir-10kbmerge-peaks.bed 

